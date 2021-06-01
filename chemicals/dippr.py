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

This module contains implementations of various numered property equations
used by the DIPPR, the Design Institude for Physical Property Research.

No actual data is included in this module; it is just functional
implementations of the formulas and some of their derivatives/integrals.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Equations
---------
.. autofunction:: chemicals.dippr.EQ100
.. autofunction:: chemicals.dippr.EQ101
.. autofunction:: chemicals.dippr.EQ102
.. autofunction:: chemicals.dippr.EQ104
.. autofunction:: chemicals.dippr.EQ105
.. autofunction:: chemicals.dippr.EQ106
.. autofunction:: chemicals.dippr.EQ107
.. autofunction:: chemicals.dippr.EQ114
.. autofunction:: chemicals.dippr.EQ115
.. autofunction:: chemicals.dippr.EQ116
.. autofunction:: chemicals.dippr.EQ127

Jacobians (for fitting)
-----------------------
.. autofunction:: chemicals.dippr.EQ101_fitting_jacobian
.. autofunction:: chemicals.dippr.EQ102_fitting_jacobian
.. autofunction:: chemicals.dippr.EQ105_fitting_jacobian
.. autofunction:: chemicals.dippr.EQ106_fitting_jacobian
.. autofunction:: chemicals.dippr.EQ107_fitting_jacobian

"""

from __future__ import division

__all__ = ['EQ100', 'EQ101', 'EQ102', 'EQ104', 'EQ105', 'EQ106', 'EQ107',
           'EQ114', 'EQ115', 'EQ116', 'EQ127', 
           'EQ101_fitting_jacobian', 'EQ102_fitting_jacobian',
           'EQ106_fitting_jacobian', 'EQ105_fitting_jacobian',
           'EQ107_fitting_jacobian']

from chemicals.utils import log, exp, sinh, cosh, atan, atanh, sqrt, tanh
from cmath import log as clog
from cmath import sqrt as csqrt
from fluids.numerics import hyp2f1, trunc_exp, trunc_log, numpy as np

order_not_found_msg = ('Only the actual property calculation, first temperature '
                       'derivative, first temperature integral, and first '
                       'temperature integral over temperature are supported '
                       'with order=  0, 1, -1, or -1j respectively')

order_not_found_pos_only_msg = ('Only the actual property calculation, and'
                                'temperature derivative(s) are supported')

def EQ100(T, A=0, B=0, C=0, D=0, E=0, F=0, G=0, order=0):
    r'''DIPPR Equation # 100. Used in calculating the molar heat capacities
    of liquids and solids, liquid thermal conductivity, and solid density.
    All parameters default to zero. As this is a straightforward polynomial,
    no restrictions on parameters apply. Note that high-order polynomials like
    this may need large numbers of decimal places to avoid unnecessary error.

    .. math::
        Y = A + BT + CT^2 + DT^3 + ET^4 + FT^5 + GT^6

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-G : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T, integral with respect to T, and integral
    over T with respect to T are computed as follows. All derivatives and
    integrals are easily computed with SymPy.

    .. math::
        \frac{d Y}{dT} = B + 2 C T + 3 D T^{2} + 4 E T^{3} + 5 F T^{4}
        + 6 G T^{5}

    .. math::
        \int Y dT = A T + \frac{B T^{2}}{2} + \frac{C T^{3}}{3} + \frac{D
        T^{4}}{4} + \frac{E T^{5}}{5} + \frac{F T^{6}}{6} + \frac{G T^{7}}{7}

    .. math::
        \int \frac{Y}{T} dT = A \ln{\left (T \right )} + B T + \frac{C T^{2}}
        {2} + \frac{D T^{3}}{3} + \frac{E T^{4}}{4} + \frac{F T^{5}}{5}
        + \frac{G T^{6}}{6}

    Examples
    --------
    Water liquid heat capacity; DIPPR coefficients normally listed in J/kmol/K.

    >>> EQ100(300, 276370., -2090.1, 8.125, -0.014116, 0.0000093701)
    75355.81000000003

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        return A + T*(B + T*(C + T*(D + T*(E + T*(F + G*T)))))
    elif order == 1:
        return B + T*(2.0*C + T*(3.0*D + T*(4.0*E + T*(5.0*F + 6.0*G*T))))
    elif order == -1:
        return T*(A + T*(B*0.5 + T*(C*(1.0/3.0) + T*(D*0.25 + T*(E*0.2 + T*(F*(1.0/6.0) + G*T*(1.0/7.0)))))))
    elif order == -1j:
        return A*log(T) + T*(B + T*(C*0.5 + T*(D*(1.0/3.0) + T*(E*0.25 + T*(F*0.2 + G*T*(1.0/6.0))))))
    else:
        raise ValueError(order_not_found_msg)


def EQ101(T, A, B, C=0.0, D=0.0, E=0.0, order=0):
    r'''DIPPR Equation # 101. Used in calculating vapor pressure, sublimation
    pressure, and liquid viscosity.
    All 5 parameters are required. E is often an integer. As the model is
    exponential, a sufficiently high temperature will cause an OverflowError.
    A negative temperature (or just low, if fit poorly) may cause a math domain
    error.

    .. math::
        Y = \exp\left(A + \frac{B}{T} + C\cdot \ln T + D \cdot T^E\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for `n`, the `nth` derivative of the property is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific]

    Notes
    -----
    This function is not integrable for either dT or Y/T dT.

    .. math::
        \frac{d Y}{dT} = \left(- \frac{B}{T^{2}} + \frac{C}{T}
        + \frac{D E T^{E}}{T}\right) e^{A + \frac{B}{T}
        + C \log{\left(T \right)} + D T^{E}}

    .. math::
        \frac{d^2 Y}{dT^2} = \frac{\left(\frac{2 B}{T} - C + D E^{2} T^{E}
        - D E T^{E} + \left(- \frac{B}{T} + C + D E T^{E}\right)^{2}\right)
        e^{A + \frac{B}{T} + C \log{\left(T \right)} + D T^{E}}}{T^{2}}

    .. math::
        \frac{d^3 Y}{dT^3} = \frac{\left(- \frac{6 B}{T} + 2 C + D E^{3} T^{E}
        - 3 D E^{2} T^{E} + 2 D E T^{E} + \left(- \frac{B}{T} + C
        + D E T^{E}\right)^{3} + 3 \left(- \frac{B}{T} + C + D E T^{E}\right)
        \left(\frac{2 B}{T} - C + D E^{2} T^{E} - D E T^{E}\right)\right)
        e^{A + \frac{B}{T} + C \log{\left(T \right)} + D T^{E}}}{T^{3}}

    Examples
    --------
    Water vapor pressure; DIPPR coefficients normally listed in Pa.

    >>> EQ101(300, 73.649, -7258.2, -7.3037, 4.1653E-6, 2)
    3537.44834545549

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    T_inv = 1.0/T
    try:
        T_E = T**E
    except:
        T_E = 1e250
    expr = trunc_exp(A + B*T_inv + C*trunc_log(T) + D*T_E)
    if order == 0:
        return expr
    elif order == 1:
        return T_inv*expr*(-B*T_inv + C + D*E*T_E)
    elif order == 2:
        x0 = (-B*T_inv + C + D*E*T_E)
        return expr*(2.0*B*T_inv - C + D*E*T_E*(E - 1.0) + x0*x0)*T_inv*T_inv
    elif order == 3:
        E2 = E*E
        E3 = E2*E
        x0 = (-B*T_inv + C + D*E*T_E)
        return expr*(-6.0*B*T_inv + 2.0*C + D*E3*T_E - 3*D*E2*T_E + 2.0*D*E*T_E
                     + x0*x0*x0
                     + 3.0*(-B*T_inv + C + D*E*T_E)*(2.0*B*T_inv - C + D*E2*T_E - D*E*T_E))*T_inv*T_inv*T_inv
    else:
        raise ValueError(order_not_found_pos_only_msg)

def EQ102(T, A, B, C=0.0, D=0.0, order=0):
    r'''DIPPR Equation # 102. Used in calculating vapor viscosity, vapor
    thermal conductivity, and sometimes solid heat capacity. High values of B
    raise an OverflowError.
    All 4 parameters are required. C and D are often 0.

    .. math::
        Y = \frac{A\cdot T^B}{1 + \frac{C}{T} + \frac{D}{T^2}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-D : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T, integral with respect to T, and integral
    over T with respect to T are computed as follows. The first derivative is
    easily computed; the two integrals required Rubi to perform the integration.

    .. math::
        \frac{d Y}{dT} = \frac{A B T^{B}}{T \left(\frac{C}{T} + \frac{D}{T^{2}}
        + 1\right)} + \frac{A T^{B} \left(\frac{C}{T^{2}} + \frac{2 D}{T^{3}}
        \right)}{\left(\frac{C}{T} + \frac{D}{T^{2}} + 1\right)^{2}}

    .. math::
        \int Y dT = - \frac{2 A T^{B + 3} \operatorname{hyp2f1}{\left (1,B + 3,
        B + 4,- \frac{2 T}{C - \sqrt{C^{2} - 4 D}} \right )}}{\left(B + 3\right)
        \left(C + \sqrt{C^{2} - 4 D}\right) \sqrt{C^{2} - 4 D}} + \frac{2 A
        T^{B + 3} \operatorname{hyp2f1}{\left (1,B + 3,B + 4,- \frac{2 T}{C
        + \sqrt{C^{2} - 4 D}} \right )}}{\left(B + 3\right) \left(C
        - \sqrt{C^{2} - 4 D}\right) \sqrt{C^{2} - 4 D}}

    .. math::
        \int \frac{Y}{T} dT = - \frac{2 A T^{B + 2} \operatorname{hyp2f1}{\left
        (1,B + 2,B + 3,- \frac{2 T}{C + \sqrt{C^{2} - 4 D}} \right )}}{\left(B
        + 2\right) \left(C + \sqrt{C^{2} - 4 D}\right) \sqrt{C^{2} - 4 D}}
        + \frac{2 A T^{B + 2} \operatorname{hyp2f1}{\left (1,B + 2,B + 3,
        - \frac{2 T}{C - \sqrt{C^{2} - 4 D}} \right )}}{\left(B + 2\right)
        \left(C - \sqrt{C^{2} - 4 D}\right) \sqrt{C^{2} - 4 D}}

    Examples
    --------
    Water vapor viscosity; DIPPR coefficients normally listed in Pa*s.

    >>> EQ102(300, 1.7096E-8, 1.1146, 0, 0)
    9.860384711890639e-06

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        easy = A/(1. + C/T + D/(T*T))
        try:
            return easy*T**B
        except:
            return 1e308
    elif order == 1:
        return (A*B*T**B/(T*(C/T + D/T**2 + 1))
                + A*T**B*(C/T**2 + 2*D/T**3)/(C/T + D/T**2 + 1)**2)
    elif order == -1: # numba: delete
        # numba-scipy does not support complex numbers so this does not work in numba
        # imaginary part is 0
        # Hours spend trying to make hyp2f1 use real inputs only: 3
        c0 = 3.0 + B # numba: delete
        x0 = csqrt(C*C - 4.0*D) # numba: delete
        arg3_hyp = (-2.0*T/(C - x0)) # numba: delete
        hyp2f1_term1 = hyp2f1(1.0, c0, 4.0+B, arg3_hyp) # numba: delete
        hyp2f1_term2 = hyp2f1_term1.real - hyp2f1_term1.imag*1.0j# numba: delete
        x5 = 2.0*A*T**(c0)/(c0*x0) # numba: delete
        x10 = x5/(C - x0) # numba: delete
        x11 = x5/(C + x0) # numba: delete
        return float((hyp2f1_term1*x10 - hyp2f1_term2*x11).real) # numba: delete
    elif order == -1j: # numba: delete
        return float((2*A*T**(2+B)*hyp2f1(1.0, 2.0+B, 3.0+B, -2*T/(C - csqrt(C*C - 4*D)))/( # numba: delete
                (2+B)*(C - csqrt(C*C-4*D))*csqrt(C*C-4*D)) -2*A*T**(2+B)*hyp2f1( # numba: delete
                1.0, 2.0+B, 3.0+B, -2*T/(C + csqrt(C*C - 4*D)))/((2+B)*(C + csqrt( # numba: delete
                C*C-4*D))*csqrt(C*C-4*D))).real) # numba: delete
    else:
        raise ValueError(order_not_found_msg)

def EQ101_fitting_jacobian(Ts, A, B, C, D, E):
    r'''Compute and return the Jacobian of the property predicted by 
    DIPPR Equation # 101 with respect to all the coefficients. This is used in
    fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]

    Returns
    -------
    jac : list[list[float, 5], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]

    '''
    N = len(Ts)
#    out = np.zeros((N, 5)) # numba: uncomment
    out = [[0.0]*5 for _ in range(N)] # numba: delete
    for i in range(N):
        x0 = log(Ts[i])
        x1 = 1.0/Ts[i]
        x2 = Ts[i]**E
        x3 = D*x2
        x4 = exp(A + B*x1 + C*x0 + x3)
        x5 = x0*x4
        out[i][0] = x4
        out[i][1] = x1*x4
        out[i][2] = x5
        out[i][3] = x2*x4
        out[i][4] = x3*x5
    return out

def EQ102_fitting_jacobian(Ts, A, B, C, D):
    r'''Compute and return the Jacobian of the property predicted by 
    DIPPR Equation # 102 with respect to all the coefficients. This is used in
    fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    A-D : float
        Parameter for the equation; chemical and property specific [-]

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
        x0 = Ts[i]**B
        x1 = 1.0/Ts[i]
        x2 = x1*x1
        x3 = C*x1 + D*x2 + 1.0
        x4 = x0/x3
        x5 = A*x4/x3
        lnT = log(Ts[i])
        out[i][0] = x4
        out[i][1] = A*x4*lnT
        out[i][2] = -x1*x5
        out[i][3] = -x2*x5
    return out

def EQ105_fitting_jacobian(Ts, A, B, C, D):
    r'''Compute and return the Jacobian of the property predicted by 
    DIPPR Equation # 105 with respect to all the coefficients. This is used in
    fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    A-D : float
        Parameter for the equation; chemical and property specific [-]

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
        r = out[i]
        x0 = 1.0 - Ts[i]/C
        
        if D < 1.0 and x0 < 0.0:
            r[0] = 1.0/B
            r[1] = -A/(B*B)
        else:
            x1 = x0**D
            x2 = x1 + 1.0
            x3 = A*B**(-x1 - 1.0)
            x4 = x1*x3*log(B)
            r[0] = B**(-x2)
            r[1] = -x2*x3/B
            r[2] = -D*Ts[i]*x4/(C*C*x0)
            r[3] = -x4*log(x0)
    return out

def EQ106_fitting_jacobian(Ts, Tc, A, B, C, D, E):
    r'''Compute and return the Jacobian of the property predicted by 
    DIPPR Equation # 106 with respect to all the coefficients. This is used in
    fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    Tc : float
        Critical temperature, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]

    Returns
    -------
    jac : list[list[float, 5], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]

    '''
    N = len(Ts)
#    out = np.zeros((N, 5)) # numba: uncomment
    out = [[0.0]*5 for _ in range(N)] # numba: delete
    for i in range(N):
        x0 = Ts[i]/Tc
        if x0 != 1.0:
            x1 = 1.0 - x0
            x2 = x1**(B + x0*(C + x0*(D + E*x0)))
            x3 = A*x2*log(x1)
            r = out[i]
            r[0] = x2
            r[1] = x3
            r[2] = x0*x3
            r[3] = x0*x0*x3
            r[4] = x0*x0*x0*x3
    return out

def EQ107_fitting_jacobian(Ts, A, B, C, D, E):
    r'''Compute and return the Jacobian of the property predicted by 
    DIPPR Equation # 107 with respect to all the coefficients. This is used in
    fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]

    Returns
    -------
    jac : list[list[float, 5], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]

    '''
    N = len(Ts)
#    out = np.zeros((N, 5)) # numba: uncomment
    out = [[0.0]*5 for _ in range(N)] # numba: delete
    for i in range(N):
        r = out[i]
        x1 = 1.0/Ts[i]
        x0 = x1*x1
        x2 = C*x1
        x3 = sinh(x2)
        x3_inv = 1.0/x3
        x4 = x0*x3_inv*x3_inv
        x5 = E*x1
        x6 = cosh(x5)
        x6_inv = 1.0/x6
        x7 = x0*x6_inv*x6_inv
        r[0] = 1.0
        r[1] = C*C*x4
        r[2] = 2.0*B*C*x4*(-x2*cosh(x2)*x3_inv + 1.0)
        r[3] = E*E*x7
        r[4] = 2.0*D*E*x7*(-x5*sinh(x5)*x6_inv + 1.0)
    return out

def EQ104(T, A, B, C=0.0, D=0.0, E=0.0, order=0):
    r'''DIPPR Equation #104. Often used in calculating second virial
    coefficients of gases. All 5 parameters are required.
    C, D, and E are normally large values.

    .. math::
        Y = A + \frac{B}{T} + \frac{C}{T^3} + \frac{D}{T^8} + \frac{E}{T^9}

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T, integral with respect to T, and integral
    over T with respect to T are computed as follows. All expressions can be
    obtained with SymPy readily.

    .. math::
        \frac{d Y}{dT} = - \frac{B}{T^{2}} - \frac{3 C}{T^{4}}
        - \frac{8 D}{T^{9}} - \frac{9 E}{T^{10}}

    .. math::
        \int Y dT = A T + B \ln{\left (T \right )} - \frac{1}{56 T^{8}}
        \left(28 C T^{6} + 8 D T + 7 E\right)

    .. math::
        \int \frac{Y}{T} dT = A \ln{\left (T \right )} - \frac{1}{72 T^{9}}
        \left(72 B T^{8} + 24 C T^{6} + 9 D T + 8 E\right)

    Examples
    --------
    Water second virial coefficient; DIPPR coefficients normally dimensionless.

    >>> EQ104(300, 0.02222, -26.38, -16750000, -3.894E19, 3.133E21)
    -1.1204179007265156

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        T2 = T*T
        return A + (B + (C + (D + E/T)/(T2*T2*T))/T2)/T
    elif order == 1:
        T2 = T*T
        T4 = T2*T2
        return (-B + (-3*C + (-8*D - 9*E/T)/(T4*T))/T2)/T2
    elif order == -1:
        return A*T + B*log(T) - (28*C*T**6 + 8*D*T + 7*E)/(56*T**8)
    elif order == -1j:
        return A*log(T) - (72*B*T**8 + 24*C*T**6 + 9*D*T + 8*E)/(72*T**9)
    else:
        raise ValueError(order_not_found_msg)


def EQ105(T, A, B, C, D, order=0):
    r'''DIPPR Equation #105. Often used in calculating liquid molar density.
    All 4 parameters are required. C is sometimes the fluid's critical
    temperature.

    .. math::
        Y = \frac{A}{B^{1 + \left(1-\frac{T}{C}\right)^D}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-D : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, 2, and 3, that derivative of the property is returned; No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific]

    Notes
    -----
    This expression can be integrated in terms of the incomplete gamma function
    for dT, however nans are the only output from that function.
    For Y/T dT no integral could be found.

    .. math::
        \frac{d Y}{dT} = \frac{A B^{- \left(1 - \frac{T}{C}\right)^{D} - 1} D
        \left(1 - \frac{T}{C}\right)^{D} \log{\left(B \right)}}{C \left(1
        - \frac{T}{C}\right)}

    .. math::
        \frac{d^2 Y}{dT^2} = \frac{A B^{- \left(1 - \frac{T}{C}\right)^{D} - 1}
        D \left(1 - \frac{T}{C}\right)^{D} \left(D \left(1 - \frac{T}{C}
        \right)^{D} \log{\left(B \right)} - D + 1\right) \log{\left(B \right)}}
        {C^{2} \left(1 - \frac{T}{C}\right)^{2}}

    .. math::
        \frac{d^3 Y}{dT^3} = \frac{A B^{- \left(1 - \frac{T}{C}\right)^{D} - 1}
        D \left(1 - \frac{T}{C}\right)^{D} \left(D^{2} \left(1 - \frac{T}{C}
        \right)^{2 D} \log{\left(B \right)}^{2} - 3 D^{2} \left(1 - \frac{T}{C}
        \right)^{D} \log{\left(B \right)} + D^{2} + 3 D \left(1 - \frac{T}{C}
        \right)^{D} \log{\left(B \right)} - 3 D + 2\right) \log{\left(B
        \right)}}{C^{3} \left(1 - \frac{T}{C}\right)^{3}}

    Examples
    --------
    Hexane molar density; DIPPR coefficients normally in kmol/m^3.

    >>> EQ105(300., 0.70824, 0.26411, 507.6, 0.27537)
    7.593170096339237

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        problematic = (1. - T/C)
        if D < 1.0 and problematic < 0.0:
            # Handle the case of a negative D exponent with a (1. - T/C) under 0 which would yield a complex number
            problematic = 0.0
        ans = A*B**(-(1. + problematic**D))
        return ans
    elif order == 1:
        x0 = 1.0/C
        x1 = 1.0 - T*x0
        x2 = x1**D
        return A*B**(-x2 - 1.0)*D*x0*x2*log(B)/x1
    elif order == 2:
        x0 = 1.0 - T/C
        x1 = x0**D
        x2 = D*x1*log(B)
        den = 1.0/(C*x0)
        return A*B**(-x1 - 1.0)*x2*(1.0 - D + x2)*den*den
    elif order == 3:
        x0 = 1.0 - T/C
        x1 = x0**D
        x2 = 3.0*D
        x3 = D*D
        x4 = log(B)
        x5 = x1*x4
        den = 1.0/(C*x0)
        return A*B**(-x1 - 1.0)*D*x5*(x0**(2.0*D)*x3*x4*x4 + x2*x5 - x2 - 3.0*x3*x5 + x3 + 2.0)*den*den*den
    else:
        raise ValueError(order_not_found_msg)



def EQ106(T, Tc, A, B, C=0.0, D=0.0, E=0.0, order=0):
    r'''DIPPR Equation #106. Often used in calculating liquid surface tension,
    and heat of vaporization.
    Only parameters A and B parameters are required; many fits include no
    further parameters. Critical temperature is also required.

    .. math::
        Y = A(1-T_r)^{B + C T_r + D T_r^2 + E T_r^3}

    .. math::
        Tr = \frac{T}{Tc}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    A-D : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, 2, and 3, that derivative of the property is returned; No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific]

    Notes
    -----
    This form is used by Yaws with only the parameters `A` and `B`.

    The integral could not be found, but the integral over T actually could,
    again in terms of hypergeometric functions.

    .. math::
        \frac{d Y}{dT} = A \left(- \frac{T}{T_{c}} + 1\right)^{B + \frac{C T}
        {T_{c}} + \frac{D T^{2}}{T_{c}^{2}} + \frac{e T^{3}}{T_{c}^{3}}} \left(
        \left(\frac{C}{T_{c}} + \frac{2 D T}{T_{c}^{2}} + \frac{3 e T^{2}}
        {T_{c}^{3}}\right) \log{\left(- \frac{T}{T_{c}} + 1 \right)} - \frac{B
        + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}} + \frac{e T^{3}}
        {T_{c}^{3}}}{T_{c} \left(- \frac{T}{T_{c}} + 1\right)}\right)

    .. math::
        \frac{d^2 Y}{dT^2} = \frac{A \left(- \frac{T}{T_{c}} + 1\right)^{B
        + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}} + \frac{e T^{3}}
        {T_{c}^{3}}} \left(2 \left(D + \frac{3 e T}{T_{c}}\right) \log{\left(
        - \frac{T}{T_{c}} + 1 \right)} + \left(\left(C + \frac{2 D T}{T_{c}}
        + \frac{3 e T^{2}}{T_{c}^{2}}\right) \log{\left(- \frac{T}{T_{c}}
        + 1 \right)} + \frac{B + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}}
        + \frac{e T^{3}}{T_{c}^{3}}}{\frac{T}{T_{c}} - 1}\right)^{2}
        + \frac{2 \left(C + \frac{2 D T}{T_{c}} + \frac{3 e T^{2}}{T_{c}^{2}}
        \right)}{\frac{T}{T_{c}} - 1} - \frac{B + \frac{C T}{T_{c}} + \frac{D
        T^{2}}{T_{c}^{2}} + \frac{e T^{3}}{T_{c}^{3}}}{\left(\frac{T}{T_{c}}
        - 1\right)^{2}}\right)}{T_{c}^{2}}

    .. math::
        \frac{d^3 Y}{dT^3} = \frac{A \left(- \frac{T}{T_{c}} + 1\right)^{B
        + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}} + \frac{e T^{3}}
        {T_{c}^{3}}} \left(\frac{6 \left(D + \frac{3 e T}{T_{c}}\right)}
        {\frac{T}{T_{c}} - 1} + \left(\left(C + \frac{2 D T}{T_{c}}
        + \frac{3 e T^{2}}{T_{c}^{2}}\right) \log{\left(- \frac{T}{T_{c}}
        + 1 \right)} + \frac{B + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}}
        + \frac{e T^{3}}{T_{c}^{3}}}{\frac{T}{T_{c}} - 1}\right)^{3}
        + 3 \left(\left(C + \frac{2 D T}{T_{c}} + \frac{3 e T^{2}}{T_{c}^{2}}
        \right) \log{\left(- \frac{T}{T_{c}} + 1 \right)} + \frac{B
        + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}} + \frac{e T^{3}}
        {T_{c}^{3}}}{\frac{T}{T_{c}} - 1}\right) \left(2 \left(D + \frac{3 e T}
        {T_{c}}\right) \log{\left(- \frac{T}{T_{c}} + 1 \right)} + \frac{2
        \left(C + \frac{2 D T}{T_{c}} + \frac{3 e T^{2}}{T_{c}^{2}}\right)}
        {\frac{T}{T_{c}} - 1} - \frac{B + \frac{C T}{T_{c}} + \frac{D T^{2}}
        {T_{c}^{2}} + \frac{e T^{3}}{T_{c}^{3}}}{\left(\frac{T}{T_{c}}
        - 1\right)^{2}}\right) + 6 e \log{\left(- \frac{T}{T_{c}} + 1 \right)}
        - \frac{3 \left(C + \frac{2 D T}{T_{c}} + \frac{3 e T^{2}}{T_{c}^{2}}
        \right)}{\left(\frac{T}{T_{c}} - 1\right)^{2}} + \frac{2 \left(B
        + \frac{C T}{T_{c}} + \frac{D T^{2}}{T_{c}^{2}} + \frac{e T^{3}}
        {T_{c}^{3}}\right)}{\left(\frac{T}{T_{c}} - 1\right)^{3}}\right)}
        {T_{c}^{3}}

    Examples
    --------
    Water surface tension; DIPPR coefficients normally in Pa*s.

    >>> EQ106(300, 647.096, 0.17766, 2.567, -3.3377, 1.9699)
    0.07231499373541

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        Tr = T/Tc
        tau = (1.0 - Tr)
        if tau == 0.0:
            return 0.0
        power = (B + Tr*(C + Tr*(D + E*Tr)))
        try:
            return A*tau**power
        except:
            # TODO: after more testing with regression, maybe return a more
            # precise value or allow A to impact the result
            return 1e300
    elif order == 1:
        x0 = 1.0/Tc
        x1 = T*x0
        x2 = 1.0 - x1
        x3 = E*x1
        x4 = C + x1*(D + x3)
        x5 = B + x1*x4
        return A*x0*x2**x5*(x5/(x1 - 1.0) + (x1*(D + 2.0*x3) + x4)*log(x2))
    elif order == 2:
        x0 = T/Tc
        x1 = 1.0 - x0
        x2 = E*x0
        x3 = C + x0*(D + x2)
        x4 = B + x0*x3
        x5 = log(x1)
        x6 = x0 - 1.0
        x7 = 1.0/x6
        x8 = x0*(D + 2.0*x2) + x3
        return (A*x1**x4*(-x4/x6**2 + 2*x5*(D + 3.0*x2) + 2.0*x7*x8
                          + (x4*x7 + x5*x8)**2)/Tc**2)
    elif order == 3:
        x0 = T/Tc
        x1 = 1.0 - x0
        x2 = E*x0
        x3 = C + x0*(D + x2)
        x4 = B + x0*x3
        x5 = log(x1)
        x6 = D + 3.0*x2
        x7 = x0 - 1.0
        x8 = 1/x7
        x9 = x7**(-2)
        x10 = x0*(D + 2.0*x2) + x3
        x11 = x10*x5 + x4*x8
        return (A*x1**x4*(-3*x10*x9 + x11**3 + 3*x11*(2*x10*x8 - x4*x9 + 2*x5*x6)
                          + 2*x4/x7**3 + 6*E*x5 + 6*x6*x8)/Tc**3)
    else:
        raise ValueError(order_not_found_msg)


def EQ107(T, A=0, B=0, C=0, D=0, E=0, order=0):
    r'''DIPPR Equation #107. Often used in calculating ideal-gas heat capacity.
    All 5 parameters are required.
    Also called the Aly-Lee equation.

    .. math::
        Y = A + B\left[\frac{C/T}{\sinh(C/T)}\right]^2 + D\left[\frac{E/T}{
        \cosh(E/T)}\right]^2

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T, integral with respect to T, and integral
    over T with respect to T are computed as follows. The derivative is
    obtained via SymPy; the integrals from Wolfram Alpha.

    .. math::
        \frac{d Y}{dT} = \frac{2 B C^{3} \cosh{\left (\frac{C}{T} \right )}}
        {T^{4} \sinh^{3}{\left (\frac{C}{T} \right )}} - \frac{2 B C^{2}}{T^{3}
        \sinh^{2}{\left (\frac{C}{T} \right )}} + \frac{2 D E^{3} \sinh{\left
        (\frac{E}{T} \right )}}{T^{4} \cosh^{3}{\left (\frac{E}{T} \right )}}
        - \frac{2 D E^{2}}{T^{3} \cosh^{2}{\left (\frac{E}{T} \right )}}

    .. math::
        \int Y dT = A T + \frac{B C}{\tanh{\left (\frac{C}{T} \right )}}
        - D E \tanh{\left (\frac{E}{T} \right )}

    .. math::
        \int \frac{Y}{T} dT = A \ln{\left (T \right )} + \frac{B C}{T \tanh{
        \left (\frac{C}{T} \right )}} - B \ln{\left (\sinh{\left (\frac{C}{T}
        \right )} \right )} - \frac{D E}{T} \tanh{\left (\frac{E}{T} \right )}
        + D \ln{\left (\cosh{\left (\frac{E}{T} \right )} \right )}

    Examples
    --------
    Water ideal gas molar heat capacity; DIPPR coefficients normally in
    J/kmol/K

    >>> EQ107(300., 33363., 26790., 2610.5, 8896., 1169.)
    33585.90452768923

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    .. [2] Aly, Fouad A., and Lloyd L. Lee. "Self-Consistent Equations for
       Calculating the Ideal Gas Heat Capacity, Enthalpy, and Entropy." Fluid
       Phase Equilibria 6, no. 3 (January 1, 1981): 169-79.
       doi:10.1016/0378-3812(81)85002-9.
    '''
    if order == 0:
        C_T = C/T
        t0 = 2.0*C_T/(trunc_exp(C_T) - trunc_exp(-C_T))
        E_T = E/T
        t1 = 2.0*E_T/(trunc_exp(-E_T) + trunc_exp(E_T))
        return A + B*t0*t0 + D*t1*t1
    elif order == 1:
        return (2*B*C**3*cosh(C/T)/(T**4*sinh(C/T)**3)
                - 2*B*C**2/(T**3*sinh(C/T)**2)
                + 2*D*E**3*sinh(E/T)/(T**4*cosh(E/T)**3)
                - 2*D*E**2/(T**3*cosh(E/T)**2))
    elif order == -1:
        return A*T + B*C/tanh(C/T) - D*E*tanh(E/T)
    elif order == -1j:
        return (A*log(T) + B*C/tanh(C/T)/T - B*log(sinh(C/T))
                - D*E*tanh(E/T)/T + D*log(cosh(E/T)))
    else:
        raise ValueError(order_not_found_msg)


def EQ114(T, Tc, A, B, C, D, order=0):
    r'''DIPPR Equation #114. Rarely used, normally as an alternate liquid
    heat capacity expression. All 4 parameters are required, as well as
    critical temperature.

    .. math::
        Y = \frac{A^2}{\tau} + B - 2AC\tau - AD\tau^2 - \frac{1}{3}C^2\tau^3
        - \frac{1}{2}CD\tau^4 - \frac{1}{5}D^2\tau^5

    .. math::
        \tau = 1 - \frac{T}{Tc}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    A-D : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T, integral with respect to T, and integral
    over T with respect to T are computed as follows. All expressions can be
    obtained with SymPy readily.

    .. math::
        \frac{d Y}{dT} = \frac{A^{2}}{T_{c} \left(- \frac{T}{T_{c}}
        + 1\right)^{2}} + \frac{2 A}{T_{c}} C + \frac{2 A}{T_{c}} D \left(
        - \frac{T}{T_{c}} + 1\right) + \frac{C^{2}}{T_{c}} \left(
        - \frac{T}{T_{c}} + 1\right)^{2} + \frac{2 C}{T_{c}} D \left(
        - \frac{T}{T_{c}} + 1\right)^{3} + \frac{D^{2}}{T_{c}} \left(
        - \frac{T}{T_{c}} + 1\right)^{4}

    .. math::
        \int Y dT = - A^{2} T_{c} \ln{\left (T - T_{c} \right )} + \frac{D^{2}
        T^{6}}{30 T_{c}^{5}} - \frac{T^{5}}{10 T_{c}^{4}} \left(C D + 2 D^{2}
        \right) + \frac{T^{4}}{12 T_{c}^{3}} \left(C^{2} + 6 C D + 6 D^{2}
        \right) - \frac{T^{3}}{3 T_{c}^{2}} \left(A D + C^{2} + 3 C D
        + 2 D^{2}\right) + \frac{T^{2}}{2 T_{c}} \left(2 A C + 2 A D + C^{2}
        + 2 C D + D^{2}\right) + T \left(- 2 A C - A D + B - \frac{C^{2}}{3}
        - \frac{C D}{2} - \frac{D^{2}}{5}\right)

    .. math::
        \int \frac{Y}{T} dT = - A^{2} \ln{\left (T + \frac{- 60 A^{2} T_{c}
        + 60 A C T_{c} + 30 A D T_{c} - 30 B T_{c} + 10 C^{2} T_{c}
        + 15 C D T_{c} + 6 D^{2} T_{c}}{60 A^{2} - 60 A C - 30 A D + 30 B
        - 10 C^{2} - 15 C D - 6 D^{2}} \right )} + \frac{D^{2} T^{5}}
        {25 T_{c}^{5}} - \frac{T^{4}}{8 T_{c}^{4}} \left(C D + 2 D^{2}
        \right) + \frac{T^{3}}{9 T_{c}^{3}} \left(C^{2} + 6 C D + 6 D^{2}
        \right) - \frac{T^{2}}{2 T_{c}^{2}} \left(A D + C^{2} + 3 C D
        + 2 D^{2}\right) + \frac{T}{T_{c}} \left(2 A C + 2 A D + C^{2}
        + 2 C D + D^{2}\right) + \frac{1}{30} \left(30 A^{2} - 60 A C
        - 30 A D + 30 B - 10 C^{2} - 15 C D - 6 D^{2}\right) \ln{\left
        (T + \frac{1}{60 A^{2} - 60 A C - 30 A D + 30 B - 10 C^{2} - 15 C D
        - 6 D^{2}} \left(- 30 A^{2} T_{c} + 60 A C T_{c} + 30 A D T_{c}
        - 30 B T_{c} + 10 C^{2} T_{c} + 15 C D T_{c} + 6 D^{2} T_{c}
        + T_{c} \left(30 A^{2} - 60 A C - 30 A D + 30 B - 10 C^{2} - 15 C D
        - 6 D^{2}\right)\right) \right )}

    Strictly speaking, the integral over T has an imaginary component, but
    only the real component is relevant and the complex part discarded.

    Examples
    --------
    Hydrogen liquid heat capacity; DIPPR coefficients normally in J/kmol/K.

    >>> EQ114(20, 33.19, 66.653, 6765.9, -123.63, 478.27)
    19423.948911676463

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        t = 1.-T/Tc
        return (A**2./t + B - 2.*A*C*t - A*D*t**2. - C**2.*t**3./3.
                - C*D*t**4./2. - D**2*t**5./5.)
    elif order == 1:
        return (A**2/(Tc*(-T/Tc + 1)**2) + 2*A*C/Tc + 2*A*D*(-T/Tc + 1)/Tc
                + C**2*(-T/Tc + 1)**2/Tc + 2*C*D*(-T/Tc + 1)**3/Tc
                + D**2*(-T/Tc + 1)**4/Tc)
    elif order == -1:
        return (-A**2*Tc*clog(T - Tc).real + D**2*T**6/(30*Tc**5)
                - T**5*(C*D + 2*D**2)/(10*Tc**4)
                + T**4*(C**2 + 6*C*D + 6*D**2)/(12*Tc**3) - T**3*(A*D + C**2
                + 3*C*D + 2*D**2)/(3*Tc**2) + T**2*(2*A*C + 2*A*D + C**2 + 2*C*D
                + D**2)/(2*Tc) + T*(-2*A*C - A*D + B - C**2/3 - C*D/2 - D**2/5))
    elif order == -1j:
        return (-A**2*clog(T + (-60*A**2*Tc + 60*A*C*Tc + 30*A*D*Tc - 30*B*Tc
                + 10*C**2*Tc + 15*C*D*Tc + 6*D**2*Tc)/(60*A**2 - 60*A*C
                - 30*A*D + 30*B - 10*C**2 - 15*C*D - 6*D**2)).real
                + D**2*T**5/(25*Tc**5) - T**4*(C*D + 2*D**2)/(8*Tc**4)
                + T**3*(C**2 + 6*C*D + 6*D**2)/(9*Tc**3) - T**2*(A*D + C**2
                + 3*C*D + 2*D**2)/(2*Tc**2) + T*(2*A*C + 2*A*D + C**2 + 2*C*D
                + D**2)/Tc + (30*A**2 - 60*A*C - 30*A*D + 30*B - 10*C**2
                - 15*C*D - 6*D**2)*clog(T + (-30*A**2*Tc + 60*A*C*Tc
                + 30*A*D*Tc - 30*B*Tc + 10*C**2*Tc + 15*C*D*Tc + 6*D**2*Tc
                + Tc*(30*A**2 - 60*A*C - 30*A*D + 30*B - 10*C**2 - 15*C*D
                - 6*D**2))/(60*A**2 - 60*A*C - 30*A*D + 30*B - 10*C**2
                - 15*C*D - 6*D**2)).real/30)
    else:
        raise ValueError(order_not_found_msg)


def EQ115(T, A, B, C=0, D=0, E=0, order=0):
    r'''DIPPR Equation #115. No major uses; has been used as an alternate
    liquid viscosity expression, and as a model for vapor pressure.
    Only parameters A and B are required.

    .. math::
        Y = \exp\left(A + \frac{B}{T} + C\ln T + D T^2 + \frac{E}{T^2}\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, 2, and 3, that derivative of the property is returned; No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific]

    Notes
    -----
    No coefficients found for this expression.
    This function is not integrable for either dT or Y/T dT.

    .. math::
        \frac{d Y}{dT} = \left(- \frac{B}{T^{2}} + \frac{C}{T} + 2 D T
        - \frac{2 E}{T^{3}}\right) e^{A + \frac{B}{T} + C \log{\left(T \right)}
        + D T^{2} + \frac{E}{T^{2}}}

    .. math::
        \frac{d^2 Y}{dT^2} = \left(\frac{2 B}{T^{3}} - \frac{C}{T^{2}} + 2 D
        + \frac{6 E}{T^{4}} + \left(\frac{B}{T^{2}} - \frac{C}{T} - 2 D T
        + \frac{2 E}{T^{3}}\right)^{2}\right) e^{A + \frac{B}{T}
        + C \log{\left(T \right)} + D T^{2} + \frac{E}{T^{2}}}

    .. math::
        \frac{d^3 Y}{dT^3} =- \left(3 \left(\frac{2 B}{T^{3}} - \frac{C}{T^{2}}
        + 2 D + \frac{6 E}{T^{4}}\right) \left(\frac{B}{T^{2}} - \frac{C}{T}
        - 2 D T + \frac{2 E}{T^{3}}\right) + \left(\frac{B}{T^{2}}
        - \frac{C}{T} - 2 D T + \frac{2 E}{T^{3}}\right)^{3} + \frac{2 \left(
        \frac{3 B}{T} - C + \frac{12 E}{T^{2}}\right)}{T^{3}}\right)
        e^{A + \frac{B}{T} + C \log{\left(T \right)} + D T^{2} + \frac{E}{T^{2}}}

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        return trunc_exp(A+B/T+C*log(T)+D*T**2 + E/T**2)
    elif order == 1:
        x0 = T**2
        x1 = 1/x0
        x2 = 1/T
        return (-(B*x1 - C*x2 - 2*D*T + 2*E/T**3)*exp(A + B*x2 + C*log(T) + D*x0 + E*x1))
    elif order == 2:
        x0 = 1/T
        x1 = T**2
        x2 = 1/x1
        x3 = 2*D
        x4 = 2/T**3
        return (B*x4 - C*x2 + 6*E/T**4 + x3 + (B*x2 - C*x0 + E*x4 - T*x3)**2)*exp(A + B*x0 + C*log(T) + D*x1 + E*x2)
    elif order == 3:
        x0 = 1/T
        x1 = B*x0
        x2 = T**2
        x3 = 1/x2
        x4 = E*x3
        x5 = 2/T**3
        x6 = 2*D
        x7 = B*x3 - C*x0 + E*x5 - T*x6
        return (-(x5*(-C + 3*x1 + 12*x4) + x7**3 + 3*x7*(B*x5 - C*x3 + 6*E/T**4
                  + x6))*exp(A + C*log(T) + D*x2 + x1 + x4))
    else:
        raise ValueError(order_not_found_msg)


def EQ116(T, Tc, A, B, C, D, E, order=0):
    r'''DIPPR Equation #116. Used to describe the molar density of water fairly
    precisely; no other uses listed. All 5 parameters are needed, as well as
    the critical temperature.

    .. math::
        Y = A + B\tau^{0.35} + C\tau^{2/3} + D\tau + E\tau^{4/3}

    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    A-E : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T and integral with respect to T are
    computed as follows. The integral divided by T with respect to T has an
    extremely complicated (but still elementary) integral which can be read
    from the source. It was computed with Rubi; the other expressions can
    readily be obtained with SymPy.

    .. math::
        \frac{d Y}{dT} = - \frac{7 B}{20 T_c \left(- \frac{T}{T_c} + 1\right)^{
        \frac{13}{20}}} - \frac{2 C}{3 T_c \sqrt[3]{- \frac{T}{T_c} + 1}}
        - \frac{D}{T_c} - \frac{4 E}{3 T_c} \sqrt[3]{- \frac{T}{T_c} + 1}

    .. math::
        \int Y dT = A T - \frac{20 B}{27} T_c \left(- \frac{T}{T_c} + 1\right)^{
        \frac{27}{20}} - \frac{3 C}{5} T_c \left(- \frac{T}{T_c} + 1\right)^{
        \frac{5}{3}} + D \left(- \frac{T^{2}}{2 T_c} + T\right) - \frac{3 E}{7}
        T_c \left(- \frac{T}{T_c} + 1\right)^{\frac{7}{3}}

    Examples
    --------
    Water liquid molar density; DIPPR coefficients normally in kmol/m^3.

    >>> EQ116(300., 647.096, 17.863, 58.606, -95.396, 213.89, -141.26)
    55.17615446406527

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        tau = 1-T/Tc
        return A + B*tau**0.35 + C*tau**(2/3.) + D*tau + E*tau**(4/3.)
    elif order == 1:
        return (-7*B/(20*Tc*(-T/Tc + 1)**(13/20))
                - 2*C/(3*Tc*(-T/Tc + 1)**(1/3))
                - D/Tc - 4*E*(-T/Tc + 1)**(1/3)/(3*Tc))
    elif order == -1:
        return (A*T - 20*B*Tc*(-T/Tc + 1)**(27/20)/27
                - 3*C*Tc*(-T/Tc + 1)**(5/3)/5 + D*(-T**2/(2*Tc) + T)
                - 3*E*Tc*(-T/Tc + 1)**(7/3)/7)
    elif order == -1j:
        # 3x increase in speed - cse via sympy
        x0 = log(T)
        x1 = 0.5*x0
        x2 = 1/Tc
        x3 = T*x2
        x4 = -x3 + 1
        x5 = 1.5*C
        x6 = x4**0.333333333333333
        x7 = 2*B
        x8 = x4**0.05
        x9 = log(-x6 + 1)
        x10 = sqrt(3)
        x11 = x10*atan(x10*(2*x6 + 1)/3)
        x12 = sqrt(5)
        x13 = 0.5*x12
        x14 = x13 + 0.5
        x15 = B*x14
        x16 = sqrt(x13 + 2.5)
        x17 = 2*x8
        x18 = -x17
        x19 = -x13
        x20 = x19 + 0.5
        x21 = B*x20
        x22 = sqrt(x19 + 2.5)
        x23 = B*x16
        x24 = 0.5*sqrt(0.1*x12 + 0.5)
        x25 = x12 + 1
        x26 = 4*x8
        x27 = -x26
        x28 = sqrt(10)*B/sqrt(x12 + 5)
        x29 = 2*x12
        x30 = sqrt(x29 + 10)
        x31 = 1/x30
        x32 = -x12 + 1
        x33 = 0.5*B*x22
        x34 = -x2*(T - Tc)
        x35 = 2*x34**0.1
        x36 = x35 + 2
        x37 = x34**0.05
        x38 = x30*x37
        x39 = 0.5*B*x16
        x40 = x37*sqrt(-x29 + 10)
        x41 = 0.25*x12
        x42 = B*(-x41 + 0.25)
        x43 = x12*x37
        x44 = x35 + x37 + 2
        x45 = B*(x41 + 0.25)
        x46 = -x43
        x47 = x35 - x37 + 2
        return A*x0 + 2.85714285714286*B*x4**0.35 - C*x1 + C*x11 + D*x0 - D*x3 - E*x1 - E*x11 + 0.75*E*x4**1.33333333333333 + 3*E*x6 + 1.5*E*x9 - x15*atan(x14*(x16 + x17)) + x15*atan(x14*(x16 + x18)) - x21*atan(x20*(x17 + x22)) + x21*atan(x20*(x18 + x22)) + x23*atan(x24*(x25 + x26)) - x23*atan(x24*(x25 + x27)) - x28*atan(x31*(x26 + x32)) + x28*atan(x31*(x27 + x32)) - x33*log(x36 - x38) + x33*log(x36 + x38) + x39*log(x36 - x40) - x39*log(x36 + x40) + x4**0.666666666666667*x5 - x42*log(x43 + x44) + x42*log(x46 + x47) + x45*log(x43 + x47) - x45*log(x44 + x46) + x5*x9 + x7*atan(x8) - x7*atanh(x8)
    else:
        raise ValueError(order_not_found_msg)


def EQ127(T, A, B, C, D, E, F, G, order=0):
    r'''DIPPR Equation #127. Rarely used, and then only in calculating
    ideal-gas heat capacity. All 7 parameters are required.

    .. math::
        Y = A+B\left[\frac{\left(\frac{C}{T}\right)^2\exp\left(\frac{C}{T}
        \right)}{\left(\exp\frac{C}{T}-1 \right)^2}\right]
        +D\left[\frac{\left(\frac{E}{T}\right)^2\exp\left(\frac{E}{T}\right)}
        {\left(\exp\frac{E}{T}-1 \right)^2}\right]
        +F\left[\frac{\left(\frac{G}{T}\right)^2\exp\left(\frac{G}{T}\right)}
        {\left(\exp\frac{G}{T}-1 \right)^2}\right]

    Parameters
    ----------
    T : float
        Temperature, [K]
    A-G : float
        Parameter for the equation; chemical and property specific [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -1j, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == -1j, unchanged from default]

    Notes
    -----
    The derivative with respect to T, integral with respect to T, and integral
    over T with respect to T are computed as follows. All expressions can be
    obtained with SymPy readily.

    .. math::
        \frac{d Y}{dT} = - \frac{B C^{3} e^{\frac{C}{T}}}{T^{4}
        \left(e^{\frac{C}{T}} - 1\right)^{2}} + \frac{2 B C^{3}
        e^{\frac{2 C}{T}}}{T^{4} \left(e^{\frac{C}{T}} - 1\right)^{3}}
        - \frac{2 B C^{2} e^{\frac{C}{T}}}{T^{3} \left(e^{\frac{C}{T}}
        - 1\right)^{2}} - \frac{D E^{3} e^{\frac{E}{T}}}{T^{4}
        \left(e^{\frac{E}{T}} - 1\right)^{2}} + \frac{2 D E^{3}
        e^{\frac{2 E}{T}}}{T^{4} \left(e^{\frac{E}{T}} - 1\right)^{3}}
        - \frac{2 D E^{2} e^{\frac{E}{T}}}{T^{3} \left(e^{\frac{E}{T}}
        - 1\right)^{2}} - \frac{F G^{3} e^{\frac{G}{T}}}{T^{4}
        \left(e^{\frac{G}{T}} - 1\right)^{2}} + \frac{2 F G^{3}
        e^{\frac{2 G}{T}}}{T^{4} \left(e^{\frac{G}{T}} - 1\right)^{3}}
        - \frac{2 F G^{2} e^{\frac{G}{T}}}{T^{3} \left(e^{\frac{G}{T}}
        - 1\right)^{2}}

    .. math::
        \int Y dT = A T + \frac{B C^{2}}{C e^{\frac{C}{T}} - C}
        + \frac{D E^{2}}{E e^{\frac{E}{T}} - E}
        + \frac{F G^{2}}{G e^{\frac{G}{T}} - G}

    .. math::
        \int \frac{Y}{T} dT = A \ln{\left (T \right )} + B C^{2} \left(
        \frac{1}{C T e^{\frac{C}{T}} - C T} + \frac{1}{C T} - \frac{1}{C^{2}}
        \ln{\left (e^{\frac{C}{T}} - 1 \right )}\right) + D E^{2} \left(
        \frac{1}{E T e^{\frac{E}{T}} - E T} + \frac{1}{E T} - \frac{1}{E^{2}}
        \ln{\left (e^{\frac{E}{T}} - 1 \right )}\right) + F G^{2} \left(
        \frac{1}{G T e^{\frac{G}{T}} - G T} + \frac{1}{G T} - \frac{1}{G^{2}}
        \ln{\left (e^{\frac{G}{T}} - 1 \right )}\right)

    Examples
    --------
    Ideal gas heat capacity of methanol; DIPPR coefficients normally in
    J/kmol/K

    >>> EQ127(20., 3.3258E4, 3.6199E4, 1.2057E3, 1.5373E7, 3.2122E3, -1.5318E7, 3.2122E3)
    33258.0

    References
    ----------
    .. [1] Design Institute for Physical Properties, 1996. DIPPR Project 801
       DIPPR/AIChE
    '''
    if order == 0:
        return (A+B*((C/T)**2*exp(C/T)/(exp(C/T) - 1)**2) +
            D*((E/T)**2*exp(E/T)/(exp(E/T)-1)**2) +
            F*((G/T)**2*exp(G/T)/(exp(G/T)-1)**2))
    elif order == 1:
        return (-B*C**3*exp(C/T)/(T**4*(exp(C/T) - 1)**2)
                + 2*B*C**3*exp(2*C/T)/(T**4*(exp(C/T) - 1)**3)
                - 2*B*C**2*exp(C/T)/(T**3*(exp(C/T) - 1)**2)
                - D*E**3*exp(E/T)/(T**4*(exp(E/T) - 1)**2)
                + 2*D*E**3*exp(2*E/T)/(T**4*(exp(E/T) - 1)**3)
                - 2*D*E**2*exp(E/T)/(T**3*(exp(E/T) - 1)**2)
                - F*G**3*exp(G/T)/(T**4*(exp(G/T) - 1)**2)
                + 2*F*G**3*exp(2*G/T)/(T**4*(exp(G/T) - 1)**3)
                - 2*F*G**2*exp(G/T)/(T**3*(exp(G/T) - 1)**2))
    elif order == -1:
        return (A*T + B*C**2/(C*exp(C/T) - C) + D*E**2/(E*exp(E/T) - E)
                + F*G**2/(G*exp(G/T) - G))
    elif order == -1j:
        return (A*log(T) + B*C**2*(1/(C*T*exp(C/T) - C*T) + 1/(C*T)
                - log(exp(C/T) - 1)/C**2) + D*E**2*(1/(E*T*exp(E/T) - E*T)
                + 1/(E*T) - log(exp(E/T) - 1)/E**2)
                + F*G**2*(1/(G*T*exp(G/T) - G*T) + 1/(G*T) - log(exp(G/T)
                - 1)/G**2))
    else:
        raise ValueError(order_not_found_msg)


dippr_eq_supported_orders = {
EQ100: (0, 1, -1, -1j),
EQ101: (0, 1, 2, 3),
EQ102: (0, 1, -1, -1j),
EQ104: (0, 1, -1, -1j),
EQ105: (0, 1, 2, 3),
EQ106: (0, 1, 2, 3),
EQ107: (0, 1, -1, -1j),
EQ114: (0, 1, -1, -1j),
EQ115: (0, 1, 2, 3),
EQ116: (0, 1, -1, -1j),
EQ127: (0, 1, -1, -1j),
}
