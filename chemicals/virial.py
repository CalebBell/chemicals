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

This module contains four estimation methods for second `B` virial coefficients,
two utility covnersions for when only `B` is considered, and two methods to
calculate `Z` from higher order virial expansions.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Utilities
-----------
.. autofunction:: chemicals.virial.B_to_Z
.. autofunction:: chemicals.virial.B_from_Z
.. autofunction:: chemicals.virial.Z_from_virial_density_form
.. autofunction:: chemicals.virial.Z_from_virial_pressure_form

Second Virial Correlations
--------------------------
.. autofunction:: chemicals.virial.BVirial_Pitzer_Curl
.. autofunction:: chemicals.virial.BVirial_Abbott
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_extended
.. autofunction:: chemicals.virial.BVirial_Xiang

Third Virial Correlations
-------------------------
.. autofunction:: chemicals.virial.CVirial_Orbey_Vera
.. autofunction:: chemicals.virial.CVirial_Liu_Xiang

Cross-Parameters
----------------
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_kij
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_Tcijs
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_Pcijs
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_omegaijs

"""
from __future__ import division

__all__ = ['BVirial_Pitzer_Curl', 'BVirial_Abbott', 'BVirial_Tsonopoulos',
           'BVirial_Tsonopoulos_extended', 'BVirial_Xiang',
           'B_to_Z', 'B_from_Z', 'Z_from_virial_density_form',
           'Z_from_virial_pressure_form', 'CVirial_Orbey_Vera', 'CVirial_Liu_Xiang',
           'Tarakad_Danner_virial_CSP_kij', 'Tarakad_Danner_virial_CSP_Tcijs',
           'Tarakad_Danner_virial_CSP_Pcijs', 'Tarakad_Danner_virial_CSP_omegaijs']

from fluids.numerics import numpy as np
from cmath import sqrt as csqrt
from chemicals.utils import log, sqrt, exp
from fluids.constants import R, R_inv


def B_to_Z(B, T, P):
    r'''Calculates the compressibility factor of a gas, given its
    second virial coefficient.

    .. math::
        Z = \frac{PV}{RT} = 1 + \frac{BP}{RT}

    Parameters
    ----------
    B : float
        Second virial coefficient, [m^3/mol]
    T : float
        Temperature, [K]
    P : float
        Pressure [Pa]

    Returns
    -------
    Z : float
        Compressibility factor, [-]

    Notes
    -----
    Other forms of the virial coefficient exist.

    Examples
    --------
    >>> B_to_Z(-0.0015, 300, 1E5)
    0.939863822478637

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return 1. + B*P/(R*T)


def B_from_Z(Z, T, P):
    r'''Calculates the second virial coefficient of a pure species, given the
    compressibility factor of the gas.

    .. math::
        B = \frac{RT(Z-1)}{P}

    Parameters
    ----------
    Z : float
        Compressibility factor, [-]
    T : float
        Temperature, [K]
    P : float
        Pressure [Pa]

    Returns
    -------
    B : float
        Second virial coefficient, [m^3/mol]

    Notes
    -----
    Other forms of the virial coefficient exist.

    Examples
    --------
    >>> B_from_Z(0.94, 300, 1E5)
    -0.0014966032712675846

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return (Z - 1.0)*R*T/P


def Z_from_virial_density_form(T, P, *args):
    r'''Calculates the compressibility factor of a gas given its temperature,
    pressure, and molar density-form virial coefficients. Any number of
    coefficients is supported.

    .. math::
        Z = \frac{PV}{RT} = 1 + \frac{B}{V} + \frac{C}{V^2} + \frac{D}{V^3}
        + \frac{E}{V^4} \dots

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    B to Z : float, optional
        Virial coefficients, [various]

    Returns
    -------
    Z : float
        Compressibility factor at T, P, and with given virial coefficients, [-]

    Notes
    -----
    For use with B or with B and C or with B and C and D, optimized equations
    are used to obtain the compressibility factor directly.
    If more coefficients are provided, uses numpy's roots function to solve
    this equation. This takes substantially longer as the solution is
    numerical.

    If no virial coefficients are given, returns 1, as per the ideal gas law.

    The units of each virial coefficient are as follows, where for B, n=1, and
    C, n=2, and so on.

    .. math::
        \left(\frac{\text{m}^3}{\text{mol}}\right)^n

    Examples
    --------
    >>> Z_from_virial_density_form(300, 122057.233762653, 1E-4, 1E-5, 1E-6, 1E-7)
    1.28434940526

    References
    ----------
    .. [1] Prausnitz, John M., Rudiger N. Lichtenthaler, and Edmundo Gomes de
       Azevedo. Molecular Thermodynamics of Fluid-Phase Equilibria. 3rd
       edition. Upper Saddle River, N.J: Prentice Hall, 1998.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering.
       Butterworth-Heinemann, 1985.
    '''
    l = len(args)
    if l == 1:
        return 1/2. + (4*args[0]*P + R*T)**0.5/(2*(R*T)**0.5)
#        return ((R*T*(4*args[0]*P + R*T))**0.5 + R*T)/(2*P)
    if l == 2:
        B, C = args[0], args[1]
        # A small imaginary part is ignored
        return (P*(-(3*B*R*T/P + R**2*T**2/P**2)/(3*(-1/2 + csqrt(3)*1j/2)*(-9*B*R**2*T**2/(2*P**2) - 27*C*R*T/(2*P) + csqrt(-4*(3*B*R*T/P + R**2*T**2/P**2)**(3+0j) + (-9*B*R**2*T**2/P**2 - 27*C*R*T/P - 2*R**3*T**3/P**3)**(2+0j))/2 - R**3*T**3/P**3)**(1/3.+0j)) - (-1/2 + csqrt(3)*1j/2)*(-9*B*R**2*T**2/(2*P**2) - 27*C*R*T/(2*P) + csqrt(-4*(3*B*R*T/P + R**2*T**2/P**2)**(3+0j) + (-9*B*R**2*T**2/P**2 - 27*C*R*T/P - 2*R**3*T**3/P**3)**(2+0j))/2 - R**3*T**3/P**3)**(1/3.+0j)/3 + R*T/(3*P))/(R*T)).real
    if l == 3:
        # Huge mess. Ideally sympy could optimize a function for quick python
        # execution. Derived with kate's text highlighting
        B, C, D = args[0], args[1], args[2]
        P2 = P**2
        RT = R*T
        BRT = B*RT
        T2 = T**2
        R2 = R**2
        RT23 = 3*R2*T2
        mCRT = -C*RT
        P2256 = 256*P2

        RT23P2256 = RT23/(P2256)
        big1 = (D*RT/P - (-BRT/P - RT23/(8*P2))**2/12 - RT*(mCRT/(4*P) - RT*(BRT/(16*P) + RT23P2256)/P)/P)
        big3 = (-BRT/P - RT23/(8*P2))
        big4 = (mCRT/P - RT*(BRT/(2*P) + R2*T2/(8*P2))/P)
        big5 = big3*(-D*RT/P + RT*(mCRT/(4*P) - RT*(BRT/(16*P) + RT23P2256)/P)/P)
        big2 = 2*big1/(3*(big3**3/216 - big5/6 + big4**2/16 + csqrt(big1**3/27 + (-big3**3/108 + big5/3 - big4**2/8)**2/4))**(1/3))
        big7 = 2*BRT/(3*P) - big2 + 2*(big3**3/216 - big5/6 + big4**2/16 + csqrt(big1**3/27 + (-big3**3/108 + big5/3 - big4**2/8)**2/4))**(1/3) + R2*T2/(4*P2)
        return (P*(((csqrt(big7)/2 + csqrt(4*BRT/(3*P) - (-2*C*RT/P - 2*RT*(BRT/(2*P) + R2*T2/(8*P2))/P)/csqrt(big7) + big2 - 2*(big3**3/216 - big5/6 + big4**2/16 + csqrt(big1**3/27 + (-big3**3/108 + big5/3 - big4**2/8)**2/4))**(1/3) + R2*T2/(2*P2))/2 + RT/(4*P))))/R/T).real

    size = l + 2
#    arr = np.ones(size, dtype=np.complex128) # numba: uncomment
    arr = [1.0]*size # numba: delete
    arr[-1] = -P/R/T
    for i in range(l):
        arr[-3-i] = args[i]
    solns = np.roots(arr)
    for rho in solns:
        if abs(rho.imag) < 1e-12 and rho.real > 0.0:
            return float(P/(R*T*rho.real))
    raise ValueError("Could not find real root")


def Z_from_virial_pressure_form(P, *args):
    r'''Calculates the compressibility factor of a gas given its pressure, and
    pressure-form virial coefficients. Any number of coefficients is supported.

    .. math::
        Z = \frac{Pv}{RT} = 1 + B'P + C'P^2 + D'P^3 + E'P^4 \dots

    Parameters
    ----------
    P : float
        Pressure, [Pa]
    B to Z : float, optional
        Pressure form Virial coefficients, [various]

    Returns
    -------
    Z : float
        Compressibility factor at P, and with given virial coefficients, [-]

    Notes
    -----
    Note that although this function does not require a temperature input, it
    is still dependent on it because the coefficients themselves normally are
    regressed in terms of temperature.

    The use of this form is less common than the density form. Its coefficients
    are normally indicated with the "'" suffix.

    If no virial coefficients are given, returns 1, as per the ideal gas law.

    The units of each virial coefficient are as follows, where for B, n=1, and
    C, n=2, and so on.

    .. math::
        \left(\frac{1}{\text{Pa}}\right)^n

    Examples
    --------
    >>> Z_from_virial_pressure_form(102919.99946855308, 4.032286555169439e-09, 1.6197059494442215e-13, 6.483855042486911e-19)
    1.00283753944

    References
    ----------
    .. [1] Prausnitz, John M., Rudiger N. Lichtenthaler, and Edmundo Gomes de
       Azevedo. Molecular Thermodynamics of Fluid-Phase Equilibria. 3rd
       edition. Upper Saddle River, N.J: Prentice Hall, 1998.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering.
       Butterworth-Heinemann, 1985.
    '''
    tot = 0.0
    fact = 1.0
    for i in range(len(args)):
        tot += args[i]*fact
        fact *= P
    return 1.0 + P*tot


### Second Virial Coefficients


def BVirial_Pitzer_Curl(T, Tc, Pc, omega, order=0):
    r'''Calculates the second virial coefficient using the model in [1]_.
    Designed for simple calculations.

    .. math::
        B_r=B^{(0)}+\omega B^{(1)}

    .. math::
        B^{(0)}=0.1445-0.33/T_r-0.1385/T_r^2-0.0121/T_r^3

    .. math::
        B^{(1)} = 0.073+0.46/T_r-0.5/T_r^2 -0.097/T_r^3 - 0.0073/T_r^8

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of B itself; for 1/2/3,
        the first/second/third derivative of B with respect to temperature; and
        for -1/-2, the first/second indefinite integral of B with respect to
        temperature. No other integrals or derivatives are implemented, and an
        exception will be raised if any other order is given.

    Returns
    -------
    B : float
        Second virial coefficient in density form or its integral/derivative if
        specified, [m^3/mol or m^3/mol/K^order]

    Notes
    -----
    Analytical models for derivatives and integrals are available for orders
    -2, -1, 1, 2, and 3, all obtained with SymPy.

    For first temperature derivative of B:

    .. math::
        \frac{d B^{(0)}}{dT} = \frac{33 Tc}{100 T^{2}} + \frac{277 Tc^{2}}{1000 T^{3}} + \frac{363 Tc^{3}}{10000 T^{4}}

    .. math::
        \frac{d B^{(1)}}{dT} = - \frac{23 Tc}{50 T^{2}} + \frac{Tc^{2}}{T^{3}} + \frac{291 Tc^{3}}{1000 T^{4}} + \frac{73 Tc^{8}}{1250 T^{9}}

    For the second temperature derivative of B:

    .. math::
        \frac{d^2 B^{(0)}}{dT^2} = - \frac{3 Tc}{5000 T^{3}} \left(1100 + \frac{1385 Tc}{T} + \frac{242 Tc^{2}}{T^{2}}\right)

    .. math::
        \frac{d^2 B^{(1)}}{dT^2} = \frac{Tc}{T^{3}} \left(\frac{23}{25} - \frac{3 Tc}{T} - \frac{291 Tc^{2}}{250 T^{2}} - \frac{657 Tc^{7}}{1250 T^{7}}\right)

    For the third temperature derivative of B:

    .. math::
        \frac{d^3 B^{(0)}}{dT^3} = \frac{3 Tc}{500 T^{4}} \left(330 + \frac{554 Tc}{T} + \frac{121 Tc^{2}}{T^{2}}\right)

    .. math::
        \frac{d^3 B^{(1)}}{dT^3} = \frac{3 Tc}{T^{4}} \left(- \frac{23}{25} + \frac{4 Tc}{T} + \frac{97 Tc^{2}}{50 T^{2}} + \frac{219 Tc^{7}}{125 T^{7}}\right)

    For the first indefinite integral of B:

    .. math::
        \int{B^{(0)}} dT = \frac{289 T}{2000} - \frac{33 Tc}{100} \ln{\left (T \right )} + \frac{1}{20000 T^{2}} \left(2770 T Tc^{2} + 121 Tc^{3}\right)

    .. math::
        \int{B^{(1)}} dT = \frac{73 T}{1000} + \frac{23 Tc}{50} \ln{\left (T \right )} + \frac{1}{70000 T^{7}} \left(35000 T^{6} Tc^{2} + 3395 T^{5} Tc^{3} + 73 Tc^{8}\right)

    For the second indefinite integral of B:

    .. math::
        \int\int B^{(0)} dT dT = \frac{289 T^{2}}{4000} - \frac{33 T}{100} Tc \ln{\left (T \right )} + \frac{33 T}{100} Tc + \frac{277 Tc^{2}}{2000} \ln{\left (T \right )} - \frac{121 Tc^{3}}{20000 T}

    .. math::
        \int\int B^{(1)} dT dT = \frac{73 T^{2}}{2000} + \frac{23 T}{50} Tc \ln{\left (T \right )} - \frac{23 T}{50} Tc + \frac{Tc^{2}}{2} \ln{\left (T \right )} - \frac{1}{420000 T^{6}} \left(20370 T^{5} Tc^{3} + 73 Tc^{8}\right)

    Examples
    --------
    Example matching that in BVirial_Abbott, for isobutane.

    >>> BVirial_Pitzer_Curl(510., 425.2, 38E5, 0.193)
    -0.00020845362479301725

    References
    ----------
    .. [1] Pitzer, Kenneth S., and R. F. Curl. "The Volumetric and
       Thermodynamic Properties of Fluids. III. Empirical Equation for the
       Second Virial Coefficient1." Journal of the American Chemical Society
       79, no. 10 (May 1, 1957): 2369-70. doi:10.1021/ja01567a007.
    '''
    Tr = T/Tc
    if order == 0:
        B0 = 0.1445 - 0.33/Tr - 0.1385/Tr**2 - 0.0121/Tr**3
        B1 = 0.073 + 0.46/Tr - 0.5/Tr**2 - 0.097/Tr**3 - 0.0073/Tr**8
    elif order == 1:
        B0 = Tc*(3300*T**2 + 2770*T*Tc + 363*Tc**2)/(10000*T**4)
        B1 = Tc*(-2300*T**7 + 5000*T**6*Tc + 1455*T**5*Tc**2 + 292*Tc**7)/(5000*T**9)
    elif order == 2:
        B0 = -3*Tc*(1100*T**2 + 1385*T*Tc + 242*Tc**2)/(5000*T**5)
        B1 = Tc*(1150*T**7 - 3750*T**6*Tc - 1455*T**5*Tc**2 - 657*Tc**7)/(1250*T**10)
    elif order == 3:
        B0 = 3*Tc*(330*T**2 + 554*T*Tc + 121*Tc**2)/(500*T**6)
        B1 = 3*Tc*(-230*T**7 + 1000*T**6*Tc + 485*T**5*Tc**2 + 438*Tc**7)/(250*T**11)
    elif order == -1:
        B0 = 289*T/2000 - 33*Tc*log(T)/100 + (2770*T*Tc**2 + 121*Tc**3)/(20000*T**2)
        B1 = 73*T/1000 + 23*Tc*log(T)/50 + (35000*T**6*Tc**2 + 3395*T**5*Tc**3 + 73*Tc**8)/(70000*T**7)
    elif order == -2:
        B0 = 289*T**2/4000 - 33*T*Tc*log(T)/100 + 33*T*Tc/100 + 277*Tc**2*log(T)/2000 - 121*Tc**3/(20000*T)
        B1 = 73*T**2/2000 + 23*T*Tc*log(T)/50 - 23*T*Tc/50 + Tc**2*log(T)/2 - (20370*T**5*Tc**3 + 73*Tc**8)/(420000*T**6)
    else:
        raise ValueError('Only orders -2, -1, 0, 1, 2 and 3 are supported.')
    Br = B0 + omega*B1
    return Br*R*Tc/Pc


def BVirial_Abbott(T, Tc, Pc, omega, order=0):
    r'''Calculates the second virial coefficient using the model in [1]_.
    Simple fit to the Lee-Kesler equation.

    .. math::
        B_r=B^{(0)}+\omega B^{(1)}

    .. math::
        B^{(0)}=0.083+\frac{0.422}{T_r^{1.6}}

    .. math::
        B^{(1)}=0.139-\frac{0.172}{T_r^{4.2}}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of B itself; for 1/2/3,
        the first/second/third derivative of B with respect to temperature; and
        for -1/-2, the first/second indefinite integral of B with respect to
        temperature. No other integrals or derivatives are implemented, and an
        exception will be raised if any other order is given.

    Returns
    -------
    B : float
        Second virial coefficient in density form or its integral/derivative if
        specified, [m^3/mol or m^3/mol/K^order]

    Notes
    -----
    Analytical models for derivatives and integrals are available for orders
    -2, -1, 1, 2, and 3, all obtained with SymPy.

    For first temperature derivative of B:

    .. math::
        \frac{d B^{(0)}}{dT} = \frac{0.6752}{T \left(\frac{T}{Tc}\right)^{1.6}}

    .. math::
        \frac{d B^{(1)}}{dT} = \frac{0.7224}{T \left(\frac{T}{Tc}\right)^{4.2}}

    For the second temperature derivative of B:

    .. math::
        \frac{d^2 B^{(0)}}{dT^2} = - \frac{1.75552}{T^{2} \left(\frac{T}{Tc}\right)^{1.6}}

    .. math::
        \frac{d^2 B^{(1)}}{dT^2} = - \frac{3.75648}{T^{2} \left(\frac{T}{Tc}\right)^{4.2}}

    For the third temperature derivative of B:

    .. math::
        \frac{d^3 B^{(0)}}{dT^3} = \frac{6.319872}{T^{3} \left(\frac{T}{Tc}\right)^{1.6}}

    .. math::
        \frac{d^3 B^{(1)}}{dT^3} = \frac{23.290176}{T^{3} \left(\frac{T}{Tc}\right)^{4.2}}

    For the first indefinite integral of B:

    .. math::
        \int{B^{(0)}} dT = 0.083 T + \frac{\frac{211}{300} Tc}{\left(\frac{T}{Tc}\right)^{0.6}}

    .. math::
        \int{B^{(1)}} dT = 0.139 T + \frac{0.05375 Tc}{\left(\frac{T}{Tc}\right)^{3.2}}

    For the second indefinite integral of B:

    .. math::
        \int\int B^{(0)} dT dT = 0.0415 T^{2} + \frac{211}{120} Tc^{2} \left(\frac{T}{Tc}\right)^{0.4}

    .. math::
        \int\int B^{(1)} dT dT = 0.0695 T^{2} - \frac{\frac{43}{1760} Tc^{2}}{\left(\frac{T}{Tc}\right)^{2.2}}

    Examples
    --------
    Example is from [1]_, p. 93, and matches the result exactly, for isobutane.

    >>> BVirial_Abbott(510., 425.2, 38E5, 0.193)
    -0.00020570185009564064

    References
    ----------
    .. [1] Smith, H. C. Van Ness Joseph M. Introduction to Chemical Engineering
       Thermodynamics 4E 1987.
    '''
    Tr = T/Tc
    if order == 0:
        B0 = 0.083 - 0.422/Tr**1.6
        B1 = 0.139 - 0.172/Tr**4.2
    elif order == 1:
        B0 = 0.6752*Tr**(-1.6)/T
        B1 = 0.7224*Tr**(-4.2)/T
    elif order == 2:
        B0 = -1.75552*Tr**(-1.6)/T**2
        B1 = -3.75648*Tr**(-4.2)/T**2
    elif order == 3:
        B0 = 6.319872*Tr**(-1.6)/T**3
        B1 = 23.290176*Tr**(-4.2)/T**3
    elif order == -1:
        B0 = 0.083*T + 211/300.*Tc*(Tr)**(-0.6)
        B1 = 0.139*T + 0.05375*Tc*Tr**(-3.2)
    elif order == -2:
        B0 = 0.0415*T**2 + 211/120.*Tc**2*Tr**0.4
        B1 = 0.0695*T**2 - 43/1760.*Tc**2*Tr**(-2.2)
    else:
        raise ValueError('Only orders -2, -1, 0, 1, 2 and 3 are supported.')
    Br = B0 + omega*B1
    return Br*R*Tc/Pc


def BVirial_Tsonopoulos(T, Tc, Pc, omega, order=0):
    r'''Calculates the second virial coefficient using the model in [1]_.

    .. math::
        B_r=B^{(0)}+\omega B^{(1)}

    .. math::
        B^{(0)}= 0.1445-0.330/T_r - 0.1385/T_r^2 - 0.0121/T_r^3 - 0.000607/T_r^8

    .. math::
        B^{(1)} = 0.0637+0.331/T_r^2-0.423/T_r^3 -0.423/T_r^3 - 0.008/T_r^8

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of B itself; for 1/2/3,
        the first/second/third derivative of B with respect to temperature; and
        for -1/-2, the first/second indefinite integral of B with respect to
        temperature. No other integrals or derivatives are implemented, and an
        exception will be raised if any other order is given.

    Returns
    -------
    B : float
        Second virial coefficient in density form or its integral/derivative if
        specified, [m^3/mol or m^3/mol/K^order]

    Notes
    -----
    A more complete expression is also available, in
    BVirial_Tsonopoulos_extended.

    Analytical models for derivatives and integrals are available for orders
    -2, -1, 1, 2, and 3, all obtained with SymPy.

    For first temperature derivative of B:

    .. math::
        \frac{d B^{(0)}}{dT} = \frac{33 Tc}{100 T^{2}} + \frac{277 Tc^{2}}{1000 T^{3}} + \frac{363 Tc^{3}}{10000 T^{4}} + \frac{607 Tc^{8}}{125000 T^{9}}

    .. math::
        \frac{d B^{(1)}}{dT} = - \frac{331 Tc^{2}}{500 T^{3}} + \frac{1269 Tc^{3}}{1000 T^{4}} + \frac{8 Tc^{8}}{125 T^{9}}

    For the second temperature derivative of B:

    .. math::
        \frac{d^2 B^{(0)}}{dT^2} = - \frac{3 Tc}{125000 T^{3}} \left(27500 + \frac{34625 Tc}{T} + \frac{6050 Tc^{2}}{T^{2}} + \frac{1821 Tc^{7}}{T^{7}}\right)

    .. math::
        \frac{d^2 B^{(1)}}{dT^2} = \frac{3 Tc^{2}}{500 T^{4}} \left(331 - \frac{846 Tc}{T} - \frac{96 Tc^{6}}{T^{6}}\right)

    For the third temperature derivative of B:

    .. math::
        \frac{d^3 B^{(0)}}{dT^3} = \frac{3 Tc}{12500 T^{4}} \left(8250 + \frac{13850 Tc}{T} + \frac{3025 Tc^{2}}{T^{2}} + \frac{1821 Tc^{7}}{T^{7}}\right)

    .. math::
        \frac{d^3 B^{(1)}}{dT^3} = \frac{3 Tc^{2}}{250 T^{5}} \left(-662 + \frac{2115 Tc}{T} + \frac{480 Tc^{6}}{T^{6}}\right)

    For the first indefinite integral of B:

    .. math::
        \int{B^{(0)}} dT = \frac{289 T}{2000} - \frac{33 Tc}{100} \ln{\left (T \right )} + \frac{1}{7000000 T^{7}} \left(969500 T^{6} Tc^{2} + 42350 T^{5} Tc^{3} + 607 Tc^{8}\right)

    .. math::
        \int{B^{(1)}} dT = \frac{637 T}{10000} - \frac{1}{70000 T^{7}} \left(23170 T^{6} Tc^{2} - 14805 T^{5} Tc^{3} - 80 Tc^{8}\right)

    For the second indefinite integral of B:

    .. math::
        \int\int B^{(0)} dT dT = \frac{289 T^{2}}{4000} - \frac{33 T}{100} Tc \ln{\left (T \right )} + \frac{33 T}{100} Tc + \frac{277 Tc^{2}}{2000} \ln{\left (T \right )} - \frac{1}{42000000 T^{6}} \left(254100 T^{5} Tc^{3} + 607 Tc^{8}\right)

    .. math::
        \int\int B^{(1)} dT dT = \frac{637 T^{2}}{20000} - \frac{331 Tc^{2}}{1000} \ln{\left (T \right )} - \frac{1}{210000 T^{6}} \left(44415 T^{5} Tc^{3} + 40 Tc^{8}\right)

    Examples
    --------
    Example matching that in BVirial_Abbott, for isobutane.

    >>> BVirial_Tsonopoulos(510., 425.2, 38E5, 0.193)
    -0.00020935295404416802

    References
    ----------
    .. [1] Tsonopoulos, Constantine. "An Empirical Correlation of Second Virial
       Coefficients." AIChE Journal 20, no. 2 (March 1, 1974): 263-72.
       doi:10.1002/aic.690200209.
    '''
    Tr = T/Tc
    if order == 0:
        B0 = 0.1445 - 0.33/Tr - 0.1385/Tr**2 - 0.0121/Tr**3 - 0.000607/Tr**8
        B1 = 0.0637 + 0.331/Tr**2 - 0.423/Tr**3 - 0.008/Tr**8
    elif order == 1:
        B0 = 33*Tc/(100*T**2) + 277*Tc**2/(1000*T**3) + 363*Tc**3/(10000*T**4) + 607*Tc**8/(125000*T**9)
        B1 = -331*Tc**2/(500*T**3) + 1269*Tc**3/(1000*T**4) + 8*Tc**8/(125*T**9)
    elif order == 2:
        B0 = -3*Tc*(27500 + 34625*Tc/T + 6050*Tc**2/T**2 + 1821*Tc**7/T**7)/(125000*T**3)
        B1 = 3*Tc**2*(331 - 846*Tc/T - 96*Tc**6/T**6)/(500*T**4)
    elif order == 3:
        B0 = 3*Tc*(8250 + 13850*Tc/T + 3025*Tc**2/T**2 + 1821*Tc**7/T**7)/(12500*T**4)
        B1 = 3*Tc**2*(-662 + 2115*Tc/T + 480*Tc**6/T**6)/(250*T**5)
    elif order == -1:
        B0 = 289*T/2000. - 33*Tc*log(T)/100. + (969500*T**6*Tc**2 + 42350*T**5*Tc**3 + 607*Tc**8)/(7000000.*T**7)
        B1 = 637*T/10000. - (23170*T**6*Tc**2 - 14805*T**5*Tc**3 - 80*Tc**8)/(70000.*T**7)
    elif order == -2:
        B0 = 289*T**2/4000. - 33*T*Tc*log(T)/100. + 33*T*Tc/100. + 277*Tc**2*log(T)/2000. - (254100*T**5*Tc**3 + 607*Tc**8)/(42000000.*T**6)
        B1 = 637*T**2/20000. - 331*Tc**2*log(T)/1000. - (44415*T**5*Tc**3 + 40*Tc**8)/(210000.*T**6)
    else:
        raise ValueError('Only orders -2, -1, 0, 1, 2 and 3 are supported.')
    Br = (B0+omega*B1)
    return Br*R*Tc/Pc


def BVirial_Tsonopoulos_extended(T, Tc, Pc, omega, a=0, b=0, species_type='',
                                 dipole=0, order=0):
    r'''Calculates the second virial coefficient using the
    comprehensive model in [1]_. See the notes for the calculation of `a` and
    `b`.

    .. math::
        \frac{BP_c}{RT_c} = B^{(0)} + \omega B^{(1)} + a B^{(2)} + b B^{(3)}

    .. math::
        B^{(0)}=0.1445-0.33/T_r-0.1385/T_r^2-0.0121/T_r^3

    .. math::
        B^{(1)} = 0.0637+0.331/T_r^2-0.423/T_r^3 -0.423/T_r^3 - 0.008/T_r^8

    .. math::
        B^{(2)} = 1/T_r^6

    .. math::
        B^{(3)} = -1/T_r^8

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]
    a : float, optional
        Fit parameter, calculated based on species_type if a is not given and
        species_type matches on of the supported chemical classes.
    b : float, optional
        Fit parameter, calculated based on species_type if a is not given and
        species_type matches on of the supported chemical classes.
    species_type : str, optional
        One of .
    dipole : float
        dipole moment, optional, [Debye]
    order : int, optional
        Order of the calculation. 0 for the calculation of B itself; for 1/2/3,
        the first/second/third derivative of B with respect to temperature; and
        for -1/-2, the first/second indefinite integral of B with respect to
        temperature. No other integrals or derivatives are implemented, and an
        exception will be raised if any other order is given.

    Returns
    -------
    B : float
        Second virial coefficient in density form or its integral/derivative if
        specified, [m^3/mol or m^3/mol/K^order]

    Notes
    -----
    Analytical models for derivatives and integrals are available for orders
    -2, -1, 1, 2, and 3, all obtained with SymPy.


    To calculate `a` or `b`, the following rules are used:

    For 'simple' or 'normal' fluids:

    .. math::
        a = 0

    .. math::
        b = 0

    For 'ketone', 'aldehyde', 'alkyl nitrile', 'ether', 'carboxylic acid',
    or 'ester' types of chemicals:

    .. math::
        a = -2.14\times 10^{-4} \mu_r - 4.308 \times 10^{-21} (\mu_r)^8

    .. math::
        b = 0

    For 'alkyl halide', 'mercaptan', 'sulfide', or 'disulfide' types of
    chemicals:

    .. math::
        a = -2.188\times 10^{-4} (\mu_r)^4 - 7.831 \times 10^{-21} (\mu_r)^8

    .. math::
        b = 0

    For 'alkanol' types of chemicals (except methanol):

    .. math::
        a = 0.0878

    .. math::
        b = 0.00908 + 0.0006957 \mu_r

    For methanol:

    .. math::
        a = 0.0878

    .. math::
        b = 0.0525

    For water:

    .. math::
        a = -0.0109

    .. math::
        b = 0

    If required, the form of dipole moment used in the calculation of some
    types of `a` and `b` values is as follows:

    .. math::
        \mu_r = 100000\frac{\mu^2(Pc/101325.0)}{Tc^2}


    For first temperature derivative of B:

    .. math::
        \frac{d B^{(0)}}{dT} = \frac{33 Tc}{100 T^{2}} + \frac{277 Tc^{2}}{1000 T^{3}} + \frac{363 Tc^{3}}{10000 T^{4}} + \frac{607 Tc^{8}}{125000 T^{9}}

    .. math::
        \frac{d B^{(1)}}{dT} = - \frac{331 Tc^{2}}{500 T^{3}} + \frac{1269 Tc^{3}}{1000 T^{4}} + \frac{8 Tc^{8}}{125 T^{9}}

    .. math::
        \frac{d B^{(2)}}{dT} = - \frac{6 Tc^{6}}{T^{7}}

    .. math::
        \frac{d B^{(3)}}{dT} = \frac{8 Tc^{8}}{T^{9}}

    For the second temperature derivative of B:

    .. math::
        \frac{d^2 B^{(0)}}{dT^2} = - \frac{3 Tc}{125000 T^{3}} \left(27500 + \frac{34625 Tc}{T} + \frac{6050 Tc^{2}}{T^{2}} + \frac{1821 Tc^{7}}{T^{7}}\right)

    .. math::
        \frac{d^2 B^{(1)}}{dT^2} = \frac{3 Tc^{2}}{500 T^{4}} \left(331 - \frac{846 Tc}{T} - \frac{96 Tc^{6}}{T^{6}}\right)

    .. math::
        \frac{d^2 B^{(2)}}{dT^2} = \frac{42 Tc^{6}}{T^{8}}

    .. math::
        \frac{d^2 B^{(3)}}{dT^2} = - \frac{72 Tc^{8}}{T^{10}}

    For the third temperature derivative of B:

    .. math::
        \frac{d^3 B^{(0)}}{dT^3} = \frac{3 Tc}{12500 T^{4}} \left(8250 + \frac{13850 Tc}{T} + \frac{3025 Tc^{2}}{T^{2}} + \frac{1821 Tc^{7}}{T^{7}}\right)

    .. math::
        \frac{d^3 B^{(1)}}{dT^3} = \frac{3 Tc^{2}}{250 T^{5}} \left(-662 + \frac{2115 Tc}{T} + \frac{480 Tc^{6}}{T^{6}}\right)

    .. math::
        \frac{d^3 B^{(2)}}{dT^3} = - \frac{336 Tc^{6}}{T^{9}}

    .. math::
        \frac{d^3 B^{(3)}}{dT^3} = \frac{720 Tc^{8}}{T^{11}}

    For the first indefinite integral of B:

    .. math::
        \int{B^{(0)}} dT = \frac{289 T}{2000} - \frac{33 Tc}{100} \ln{\left (T \right )} + \frac{1}{7000000 T^{7}} \left(969500 T^{6} Tc^{2} + 42350 T^{5} Tc^{3} + 607 Tc^{8}\right)

    .. math::
        \int{B^{(1)}} dT = \frac{637 T}{10000} - \frac{1}{70000 T^{7}} \left(23170 T^{6} Tc^{2} - 14805 T^{5} Tc^{3} - 80 Tc^{8}\right)

    .. math::
        \int{B^{(2)}} dT = - \frac{Tc^{6}}{5 T^{5}}

    .. math::
        \int{B^{(3)}} dT = \frac{Tc^{8}}{7 T^{7}}

    For the second indefinite integral of B:

    .. math::
        \int\int B^{(0)} dT dT = \frac{289 T^{2}}{4000} - \frac{33 T}{100} Tc \ln{\left (T \right )} + \frac{33 T}{100} Tc + \frac{277 Tc^{2}}{2000} \ln{\left (T \right )} - \frac{1}{42000000 T^{6}} \left(254100 T^{5} Tc^{3} + 607 Tc^{8}\right)

    .. math::
        \int\int B^{(1)} dT dT = \frac{637 T^{2}}{20000} - \frac{331 Tc^{2}}{1000} \ln{\left (T \right )} - \frac{1}{210000 T^{6}} \left(44415 T^{5} Tc^{3} + 40 Tc^{8}\right)

    .. math::
        \int\int B^{(2)} dT dT = \frac{Tc^{6}}{20 T^{4}}

    .. math::
        \int\int B^{(3)} dT dT = - \frac{Tc^{8}}{42 T^{6}}

    Examples
    --------
    Example from Perry's Handbook, 8E, p2-499. Matches to a decimal place.

    >>> BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608, a=0, b=0, species_type='ketone', dipole=1.469)
    -9.679718337596426e-05

    References
    ----------
    .. [1] Tsonopoulos, C., and J. L. Heidman. "From the Virial to the Cubic
       Equation of State." Fluid Phase Equilibria 57, no. 3 (1990): 261-76.
       doi:10.1016/0378-3812(90)85126-U
    .. [2] Tsonopoulos, Constantine, and John H. Dymond. "Second Virial
       Coefficients of Normal Alkanes, Linear 1-Alkanols (and Water), Alkyl
       Ethers, and Their Mixtures." Fluid Phase Equilibria, International
       Workshop on Vapour-Liquid Equilibria and Related Properties in Binary
       and Ternary Mixtures of Ethers, Alkanes and Alkanols, 133, no. 1-2
       (June 1997): 11-34. doi:10.1016/S0378-3812(97)00058-7.
    '''
    Tr = T/Tc
    if order == 0:
        B0 = 0.1445 - 0.33/Tr - 0.1385/Tr**2 - 0.0121/Tr**3 - 0.000607/Tr**8
        B1 = 0.0637 + 0.331/Tr**2 - 0.423/Tr**3 - 0.008/Tr**8
        B2 = 1./Tr**6
        B3 = -1./Tr**8
    elif order == 1:
        B0 = 33*Tc/(100*T**2) + 277*Tc**2/(1000*T**3) + 363*Tc**3/(10000*T**4) + 607*Tc**8/(125000*T**9)
        B1 = -331*Tc**2/(500*T**3) + 1269*Tc**3/(1000*T**4) + 8*Tc**8/(125*T**9)
        B2 = -6.0*Tc**6/T**7
        B3 = 8.0*Tc**8/T**9
    elif order == 2:
        B0 = -3*Tc*(27500 + 34625*Tc/T + 6050*Tc**2/T**2 + 1821*Tc**7/T**7)/(125000*T**3)
        B1 = 3*Tc**2*(331 - 846*Tc/T - 96*Tc**6/T**6)/(500*T**4)
        B2 = 42.0*Tc**6/T**8
        B3 = -72.0*Tc**8/T**10
    elif order == 3:
        B0 = 3*Tc*(8250 + 13850*Tc/T + 3025*Tc**2/T**2 + 1821*Tc**7/T**7)/(12500*T**4)
        B1 = 3*Tc**2*(-662 + 2115*Tc/T + 480*Tc**6/T**6)/(250*T**5)
        B2 = -336.0*Tc**6/T**9
        B3 = 720.0*Tc**8/T**11
    elif order == -1:
        B0 = 289*T/2000. - 33*Tc*log(T)/100. + (969500*T**6*Tc**2 + 42350*T**5*Tc**3 + 607*Tc**8)/(7000000.*T**7)
        B1 = 637*T/10000. - (23170*T**6*Tc**2 - 14805*T**5*Tc**3 - 80*Tc**8)/(70000.*T**7)
        B2 = -Tc**6/(5*T**5)
        B3 = Tc**8/(7*T**7)
    elif order == -2:
        B0 = 289*T**2/4000. - 33*T*Tc*log(T)/100. + 33*T*Tc/100. + 277*Tc**2*log(T)/2000. - (254100*T**5*Tc**3 + 607*Tc**8)/(42000000.*T**6)
        B1 = 637*T**2/20000. - 331*Tc**2*log(T)/1000. - (44415*T**5*Tc**3 + 40*Tc**8)/(210000.*T**6)
        B2 = Tc**6/(20*T**4)
        B3 = -Tc**8/(42*T**6)
    else:
        raise ValueError('Only orders -2, -1, 0, 1, 2 and 3 are supported.')
    if a == 0 and b == 0 and species_type != '':
        if species_type == 'simple' or species_type == 'normal':
            a, b = 0, 0
        elif species_type == 'methyl alcohol':
            a, b = 0.0878, 0.0525
        elif species_type == 'water':
            a, b = -0.0109, 0
        elif dipole != 0 and Tc != 0 and Pc != 0:
            dipole_r = 1E5*dipole**2*(Pc/101325.0)/Tc**2

            if (species_type == 'ketone' or species_type == 'aldehyde'
            or species_type == 'alkyl nitrile' or species_type == 'ether'
            or species_type == 'carboxylic acid' or species_type == 'ester'):
                a, b = -2.14E-4*dipole_r-4.308E-21*dipole_r**8, 0
            elif (species_type == 'alkyl halide' or species_type == 'mercaptan'
            or species_type == 'sulfide' or species_type == 'disulfide'):
                a, b = -2.188E-4*dipole_r**4-7.831E-21*dipole_r**8, 0

            elif species_type == 'alkanol':
                a, b = 0.0878, 0.00908+0.0006957*dipole_r
    Br = B0 + omega*B1 + a*B2 + b*B3
    return Br*R*Tc/Pc

def BVirial_Xiang(T, Tc, Pc, Vc, omega):
    r'''Calculates the second virial coefficient using the model in [1]_.

    .. math::
        B = \frac{\left(-b_0T_r^{-3/4}\exp(b_1T_r^{-3}) + b_2T_r^{-1/2})
            \right)}V_c
    
    .. math::
        b_0 = b_{00} + b_{01}\omega + b_{02}\theta
        
    .. math::
        b_1 = b_{10} + b_{11}\omega + b_{12}\theta

    .. math::
        b_2 = b_{20} + b_{21}\omega + b_{22}\theta

    .. math::
        \theta = (Z_c - 0.29)^2

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    Vc : float
        Critical volume of the fluid [m^3/mol]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    B : float
        Second virial coefficient in density form [m^3/mol]
    dB_dT : float
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2 : float
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3 : float
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----

    Examples
    --------

    >>> BVirial_Xiang(388.26, 647.1, 22050000.0, 5.543076e-05, 0.344)
    (-0.0004799570, 4.6778266e-06, -7.0157656e-08, 1.4137862e-09)

    References
    ----------
    .. [1] Xiang, H. W. "The New Simple Extended Corresponding-States 
       Principle: Vapor Pressure and Second Virial Coefficient." Chemical 
       Engineering Science 57, no. 8 (April 2002): 1439049. 
       https://doi.org/10.1016/S0009-2509(02)00017-9.
    '''
    b00 = 4.553
    b01 = 4.172 
    b02 = 0.0
    b10 = 0.02644
    b11 = 0.075
    b12 = 16.5
    b20 = 3.530
    b21 = 4.297
    b22 = 0.0
    
    # TODO optimize the powers and divides
    x0 = 1.0/Tc
    x1 = T*x0
    x2 = 10000.0*omega
    x3 = (100.0*Pc*Vc*x0/R - 29.0)
    x3 = x3*x3

    x1_sqrt_inv = 1.0/sqrt(x1)
    x4 = (10000.0*b20 + b21*x2 + b22*x3)*x1_sqrt_inv
    x5 = b11*omega
    x6 = b12*x3
    T_inv = 1.0/T
    T_inv3 = T_inv*T_inv*T_inv
    Tc3 = Tc*Tc*Tc
    x8 = Tc3*T_inv3
    x9 = (10000.0*b00 + b01*x2 + b02*x3)*exp(x8*(b10 + x5 + x6*(1/10000)))*(x1_sqrt_inv*sqrt(x1_sqrt_inv))
    x10 = 3.0*x9
    x11 = 10000.0*b10 + 10000.0*x5 + x6
    x12 = x11*x8
    x13 = x12*x9
    x14 = Tc3*Tc3*x11*x11*T_inv3*T_inv3
    B = Vc*(x4 - x9)*(1/10000)
    dB = -Vc*(-x10*x12 + 5000.0*x4 - 7500.0*x9)*T_inv*(1/(100000000))
    d2B = -3.0*Vc*(x10*x14 + 55000.0*x13 - 25000000.0*x4 + 43750000.0*x9)*T_inv*T_inv*(1/(1000000000000))
    d3B = 3.0*Vc*T_inv3*(3293750000.0*x13 + 427500.*x14*x9 - 625000000000.*x4 + 1203125000000.*x9 + 9.*Tc3*Tc3*Tc3*x11*x11*x11*x9*T_inv3*T_inv3*T_inv3)*(1/10000000000000000)
    return (B, dB, d2B, d3B)

def CVirial_Orbey_Vera(T, Tc, Pc, omega):
    r'''Calculates the third virial coefficient using the model in [1]_.
    
    .. math::
        C = (RT_c/P_c)^2 (fC_{Tr}^{(0)} + \omega fC_{Tr}^{(1)})

    .. math::
        fC_{Tr}^{(0)} = 0.01407 + 0.02432T_r^{-2.8} - 0.00313T_r^{-10.5}
        
    .. math::
        fC_{Tr}^{(1)} = -0.02676 + 0.01770T_r^{-2.8} + 0.040T_r^{-3} - 0.003T_r^{-6} - 0.00228T_r^{-10.5}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    C : float
        Second virial coefficient in density form [m^6/mol^2]
    dC_dT : float
        First temperature derivative of second virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2 : float
        Second temperature derivative of second virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3 : float
        Third temperature derivative of second virial coefficient in density
        form [m^6/mol^2/K^3]
    
    Notes
    -----

    Examples
    --------
    n-octane
    
    >>> CVirial_Orbey_Vera(T=300, Tc=568.7, Pc=2490000.0, omega=0.394)
    (-1.1107124e-05, 4.1326808e-07, -1.6041435e-08, 6.7035158e-10)

    References
    ----------
    .. [1] Orbey, Hasan, and J. H. Vera. "Correlation for the Third Virial 
       Coefficient Using Tc, Pc and ω as Parameters." AIChE Journal 29, no. 1 
       (January 1, 1983): 107-13. https://doi.org/10.1002/aic.690290115.
    '''
    x0 = T/Tc
    Tinv = 1.0/T
    Tinv2 = Tinv*Tinv
    x7 = R*Tc/Pc
    x7 *= x7
    Tc3 = Tc*Tc*Tc
    x3 = Tc3*Tc3*Tinv2*Tinv2*Tinv2
    x4 = Tinv2*Tinv
    x5 = Tc3*x4
    x1 = x0**(-21.0/2.0)
    x2 = x0**(-14.0/5.0)
    x6 = -2000.0*x5
    x8 = 60.0*omega
    
    C = -x7*(2.0*omega*(114.0*x1 - 885.0*x2 + 150.0*x3 + x6 + 1338.0) + 313.0*x1 - 2432.0*x2 - 1407.0)*(1/100000)
    dC = x7*(32865.0*x1 - 68096.0*x2 + x8*(399.0*x1 - 826.0*x2 + 300.0*x3 + x6))*Tinv*(1/1000000)
    d2C = -x7*(3779475.0*x1 - 2587648.0*x2 + x8*(45885.0*x1 - 31388.0*x2 + 21000.0*x3 - 80000.0*x5))*Tinv2*(1/10000000)
    d3C = 3.0*x4*x7*(20.0*omega*(5735625.0*x1 - 1506624.0*x2 + 1680000.0*x3 - 4000000.0*x5) + 157478125.0*x1 - 41402368.0*x2)*(1/100000000)
    return C, dC, d2C, d3C

def CVirial_Liu_Xiang(T, Tc, Pc, Vc, omega):
    r'''Calculates the third virial coefficient using the model in [1]_.
    
    .. math::
        C = V_c^2 (f_{T_r}^{(0)} + \omega f_{T_r}^{(1)} + \theta f_{T_r}^{(2)})
        

    .. math::
        f_{T_r}^{(0)} = a_{00} + a_{10}T_r^{-3} + a_{20}T_r^{-6} + a_{30}T_r^{-11}
        
    .. math::
        f_{T_r}^{(1)} = a_{01} + a_{11}T_r^{-3} + a_{21}T_r^{-6} + a_{31}T_r^{-11}

    .. math::
        f_{T_r}^{(2)} = a_{02} + a_{12}T_r^{-3} + a_{22}T_r^{-6} + a_{32}T_r^{-11}
        
    .. math::
        \theta = (Z_c - 0.29)^2

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    Vc : float
        Critical volume of the fluid [m^3/mol]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    C : float
        Second virial coefficient in density form [m^6/mol^2]
    dC_dT : float
        First temperature derivative of second virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2 : float
        Second temperature derivative of second virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3 : float
        Third temperature derivative of second virial coefficient in density
        form [m^6/mol^2/K^3]
    
    Notes
    -----

    Examples
    --------
    Water at Tr = 0.6
    
    >>> CVirial_Liu_Xiang(388.26, 647.1, 22050000.0, 5.543076923076923e-05, 0.344)
    (-1.4779977e-07, 4.9949901e-09, -1.652899e-10, 5.720067e-12)

    References
    ----------
    .. [1] Liu, D. X., and H. W. Xiang. "Corresponding-States Correlation and
       Prediction of Third Virial Coefficients for a Wide Range of Substances."
       International Journal of Thermophysics 24, no. 6 (November 1, 2003):
       1667-80. https://doi.org/10.1023/B:IJOT.0000004098.98614.38.
    '''
    a00 = 0.1623538
    a01 = -0.5390344
    a02 = 34.22804
    a10 = 0.3087440
    a11 = 1.783526
    a12 = -74.76559
    a20 = -0.01790184
    a21 = -1.055391
    a22 = 279.9220
    a30 = -0.02789157
    a31 = 0.09955867
    a32 = -62.85431
    x0 = Vc*Vc
    T_inv = 1.0/T
    T_inv2 = T_inv*T_inv
    T_inv3 = T_inv2*T_inv
    T_inv8 = T_inv3*T_inv3*T_inv*T_inv
    Tc2 = Tc*Tc
    Tc3 = Tc2*Tc
    Tc8 = Tc3*Tc3*Tc2
    x2 = T_inv3*T_inv3
    x1 = Tc8*Tc3*T_inv8*T_inv3
    x3 = Tc3*Tc3*x2
    x5 = Tc3*T_inv3

    Zc = Pc*Vc*R_inv/Tc
    theta = (Zc - 0.29)
    theta *= theta
    x7 = Tc8*T_inv8
    x8 = 11.0*x7
    x9 = 6.0*x5
    x10 = x0*Tc3
    x11 = 22.0*x7
    x12 = 7.0*x5
    x13 = 143.0*x7
    x14 = 28.0*x5
    C = x0*(a00 + a10*x5 + a20*x3 + a30*x1 + omega*(a01 + a11*x5 + a21*x3 + a31*x1) + theta*(a02 + a12*x5 + a22*x3 + a32*x1))
    dC = -x10*(3.0*a10 + a20*x9 + a30*x8 + omega*(3.0*a11 + a21*x9 + a31*x8) + theta*(3.0*a12 + a22*x9 + a32*x8))*T_inv2*T_inv2
    d2C = 6.0*x10*(2.0*a10 + a20*x12 + a30*x11 + omega*(2.0*a11 + a21*x12 + a31*x11) + theta*(2.0*a12 + a22*x12 + a32*x11))*T_inv2*T_inv3
    d3C = -12.0*x10*x2*(5.0*a10 + a20*x14 + a30*x13 + omega*(5.0*a11 + a21*x14 + a31*x13) + theta*(5.0*a12 + a22*x14 + a32*x13))
    return C, dC, d2C, d3C


### Mixing Rules

def Tarakad_Danner_virial_CSP_kij(Vcs):
    r'''Calculates a binary interaction parameter for the calculation of Bij
    binary virial coefficient as shown in [1]_ and [2]_.
    
    This equation for kij is:
        
    .. math::
        k_{ij} = 1 - \frac{8\sqrt{v_{ci}v_{cj}}}{(V_{ci}^{1/3} +V_{ci}^{1/3})^3}
        
    The equation this kij is used in is 

    .. math::
        T_{cij} = \sqrt{T_{ci}T_{cj}}(1-k_{ij})

    Parameters
    ----------
    Vcs : list[float]
        Critical volumes for each species, [m^3/mol]

    Returns
    -------
    kijs : list[list[float]]
        Binary interaction parameters, [-]
    
    Notes
    -----

    Examples
    --------
    >>> Tarakad_Danner_virial_CSP_kij(Vcs=[0.000168, 0.000316])
    [[0.0, 0.01646332091], [0.0164633209, 0.0]]

    References
    ----------
    .. [1] Tarakad, Ramanathan R., and Ronald P. Danner. "An Improved 
       Corresponding States Method for Polar Fluids: Correlation of Second 
       Virial Coefficients." AIChE Journal 23, no. 5 (1977): 685-95. 
       https://doi.org/10.1002/aic.690230510.
    .. [2] Meng, Long, and Yuan-Yuan Duan. "Prediction of the Second Cross 
       Virial Coefficients of Nonpolar Binary Mixtures." Fluid Phase Equilibria
       238 (December 1, 2005): 229-38.
       https://doi.org/10.1016/j.fluid.2005.10.007.
    '''
    N = len(Vcs)
    kijs = [[0.0]*N for i in range(N)] # numba: delete
#     kijs = zeros((N, N)) # numba: uncomment
    Vc_cbrts = [0.0]*N
    for i in range(N):
        Vc_cbrts[i] = Vcs[i]**(1.0/3.0)
        
    rt8 = 2.8284271247461903 #sqrt(8)
    Vc_sqrts = [0.0]*N
    for i in range(N):
        Vc_sqrts[i] = rt8*sqrt(Vcs[i])
        
    # There is Symmetry here but it is not used
    for i in range(N):
        r = kijs[i]
        Vci_cbrt = Vc_cbrts[i]
        Vci_sqrt = Vc_sqrts[i]
        for j in range(N):
            den = Vci_cbrt + Vc_cbrts[j]
            r[j] = 1.0 - Vci_sqrt*Vc_sqrts[j]/(den*den*den)
        # More efficient and numerical error makes this non-zero
        r[i] = 0.0
    return kijs

def Tarakad_Danner_virial_CSP_Tcijs(Tcs, kijs):
    r'''Calculates the corresponding states critical temperature for the 
    calculation of Bij
    binary virial coefficient as shown in [1]_ and [2]_.
    
    .. math::
        T_{cij} = \sqrt{T_{ci}T_{cj}}(1-k_{ij})

    Parameters
    ----------
    Tcs : list[float]
        Critical temperatures for each species, [K]
    kijs : list[list[float]]
        Binary interaction parameters, [-]

    Returns
    -------
    Tcijs : list[list[float]]
        CSP Critical temperatures for each pair of species, [K]

    Notes
    -----

    Examples
    --------
    >>> kijs = Tarakad_Danner_virial_CSP_kij(Vcs=[0.000168, 0.000316])
    >>> Tarakad_Danner_virial_CSP_Tcijs(Tcs=[514.0, 591.75], kijs=kijs)
    [[514.0, 542.42694], [542.42694, 591.75000]]

    References
    ----------
    .. [1] Tarakad, Ramanathan R., and Ronald P. Danner. "An Improved 
       Corresponding States Method for Polar Fluids: Correlation of Second 
       Virial Coefficients." AIChE Journal 23, no. 5 (1977): 685-95. 
       https://doi.org/10.1002/aic.690230510.
    .. [2] Meng, Long, and Yuan-Yuan Duan. "Prediction of the Second Cross 
       Virial Coefficients of Nonpolar Binary Mixtures." Fluid Phase Equilibria
       238 (December 1, 2005): 229-38.
       https://doi.org/10.1016/j.fluid.2005.10.007.
    '''
    N = len(Tcs)
    Tc_sqrts = [0.0]*N
    for i in range(N):
        Tc_sqrts[i] = sqrt(Tcs[i])
    Tcijs = [[0.0]*N for i in range(N)] # numba: delete
#     Tcijs = zeros((N, N)) # numba: uncomment
    for i in range(N):
        # also symmetric
        kij_row = kijs[i]
        Tcij_row = Tcijs[i]
        Tci = Tc_sqrts[i]
        for j in range(N):
            Tcij_row[j] = Tci*Tc_sqrts[j]*(1.0 - kij_row[j])
    return Tcijs

def Tarakad_Danner_virial_CSP_Pcijs(Tcs, Pcs, Vcs, Tcijs):
    r'''Calculates the corresponding states critical pressure for the 
    calculation of Bij
    binary virial coefficient as shown in [1]_ and [2]_.
    
    .. math::
        P_{cij} = \frac{4T_{cij} \left(  
            \frac{P_{ci}V_{ci}}{T_{ci}} + \frac{P_{cj}V_{cj}}{T_{cj}}
            \right)
            }{(V_{ci}^{1/3} +V_{ci}^{1/3})^3}

    Parameters
    ----------
    Tcs : list[float]
        Critical temperatures for each species, [K]
    Pcs : list[float]
        Critical pressures for each species, [Pa]
    Vcs : list[float]
        Critical volumes for each species, [m^3/mol]
    Tcijs : list[list[float]]
        CSP Critical temperatures for each pair of species, [K]

    Returns
    -------
    Pcijs : list[list[float]]
        CSP Critical pressures for each pair of species, [Pa]

    Notes
    -----

    Examples
    --------
    >>> kijs = Tarakad_Danner_virial_CSP_kij(Vcs=[0.000168, 0.000316])
    >>> Tcijs = Tarakad_Danner_virial_CSP_Tcijs(Tcs=[514.0, 591.75], kijs=kijs)
    >>> Tarakad_Danner_virial_CSP_Pcijs(Tcs=[514.0, 591.75], Pcs=[6137000.0, 4108000.0], Vcs=[0.000168, 0.000316], Tcijs=Tcijs)
    [[6136999.9, 4861936.4], [4861936.4, 4107999.9]]
    
    References
    ----------
    .. [1] Tarakad, Ramanathan R., and Ronald P. Danner. "An Improved 
       Corresponding States Method for Polar Fluids: Correlation of Second 
       Virial Coefficients." AIChE Journal 23, no. 5 (1977): 685-95. 
       https://doi.org/10.1002/aic.690230510.
    .. [2] Meng, Long, and Yuan-Yuan Duan. "Prediction of the Second Cross 
       Virial Coefficients of Nonpolar Binary Mixtures." Fluid Phase Equilibria
       238 (December 1, 2005): 229-38.
       https://doi.org/10.1016/j.fluid.2005.10.007.
    '''
    N = len(Vcs)
    Pcijs = [[0.0]*N for i in range(N)] # numba: delete
#     Pcijs = zeros((N, N)) # numba: uncomment
    Vc_cbrts = [0.0]*N
    for i in range(N):
        Vc_cbrts[i] = Vcs[i]**(1.0/3.0)
    factors = [0.0]*N
    for i in range(N):
        factors[i] = 4.0*Pcs[i]*Vcs[i]/Tcs[i]
        
    for i in range(N):
        Vci_cbrt = Vc_cbrts[i]
        factori = factors[i]
        Tcij_row = Tcijs[i]
        Pcij_row = Pcijs[i]
        for j in range(N):
            den = Vci_cbrt + Vc_cbrts[j]
            Pcij_row[j] = Tcij_row[j]*(factori + factors[j])/(den*den*den)
            
#             Pcijs[i][j] = 4.0*Tcijs[i][j]*(Pcs[i]*Vcs[i]/Tcs[i]
#             + Pcs[j]*Vcs[j]/Tcs[j])/(Vcs[i]**(1/3) + Vcs[j]**(1/3)  )**3
            
    return Pcijs

def Tarakad_Danner_virial_CSP_omegaijs(omegas):
    r'''Calculates the corresponding states acentric factor for the 
    calculation of Bij
    binary virial coefficient as shown in [1]_ and [2]_.
    
    .. math::
        \omega_{ij} = 0.5(\omega_i + \omega_j)

    Parameters
    ----------
    omegas : list[float]
        Acentric factor for each species, [-]

    Returns
    -------
    omegaijs : list[list[float]]
        CSP acentric factors for each pair of species, [-]

    Notes
    -----

    Examples
    --------
    >>> Tarakad_Danner_virial_CSP_omegaijs([0.635, 0.257])
    [[0.635, 0.446], [0.446, 0.257]]
    
    References
    ----------
    .. [1] Tarakad, Ramanathan R., and Ronald P. Danner. "An Improved 
       Corresponding States Method for Polar Fluids: Correlation of Second 
       Virial Coefficients." AIChE Journal 23, no. 5 (1977): 685-95. 
       https://doi.org/10.1002/aic.690230510.
    .. [2] Meng, Long, and Yuan-Yuan Duan. "Prediction of the Second Cross 
       Virial Coefficients of Nonpolar Binary Mixtures." Fluid Phase Equilibria
       238 (December 1, 2005): 229-38.
       https://doi.org/10.1016/j.fluid.2005.10.007.
    '''
    N = len(omegas)
    omegaijs = [[0.0]*N for i in range(N)] # numba: delete
#     omegaijs = zeros((N, N)) # numba: uncomment
    for i in range(N):
        omegai = omegas[i]
        r = omegaijs[i]
        for j in range(N):
            r[j] = 0.5*(omegai + omegas[j])
    return omegaijs
