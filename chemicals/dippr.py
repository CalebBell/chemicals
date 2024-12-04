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


__all__ = ['EQ100', 'EQ101', 'EQ102', 'EQ104', 'EQ105', 'EQ106', 'EQ107',
           'EQ114', 'EQ115', 'EQ116', 'EQ127',
           'EQ101_fitting_jacobian', 'EQ102_fitting_jacobian',
           'EQ106_fitting_jacobian', 'EQ105_fitting_jacobian',
           'EQ107_fitting_jacobian',
           'EQ106_AB', 'EQ106_ABC']

from cmath import log as clog
from cmath import sqrt as csqrt
from math import atan, atanh, cosh, sinh, tanh

from fluids.numerics import exp, hyp2f1, log, sqrt, trunc_exp, trunc_log, cbrt

order_not_found_msg = ('Only the actual property calculation, first temperature '
                       'derivative, first temperature integral, and first '
                       'temperature integral over temperature are supported '
                       'with order=  0, 1, -1, or -10 respectively')

order_not_found_pos_only_msg = ('Only the actual property calculation, and'
                                'temperature derivative(s) are supported')

# Form of an enum
BASE_CALTULATION = 0
DERIVATIVE_CALCULATION = 1
SECOND_DERIVATIVE_CALCULATION = 2
THIRD_DERIVATIVE_CALCULATION = 3
FOURTH_DERIVATIVE_CALCULATION = 3
INTEGRAL_CALCULATION = -1
INTEGRAL_OVER_T_CALCULATION = -10


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
    A : float, optional
        Zero-order coefficient, default=0 [-]
    B : float, optional
        First-order coefficient, default=0 [1/K]
    C : float, optional
        Second-order coefficient, default=0 [1/K^2]
    D : float, optional
        Third-order coefficient, default=0 [1/K^3]
    E : float, optional
        Fourth-order coefficient, default=0 [1/K^4]
    F : float, optional
        Fifth-order coefficient, default=0 [1/K^5]
    G : float, optional
        Sixth-order coefficient, default=0 [1/K^6]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -10, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == INTEGRAL_OVER_T_CALCULATION, unchanged from default]

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
    elif order == INTEGRAL_OVER_T_CALCULATION:
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
    A : float
        First coefficient [-]
    B : float  
        Second coefficient [K]
    C : float, optional
        Third coefficient, default=0 [-]
    D : float, optional
        Fourth coefficient, default=0 [-]
    E : float, optional
        Fifth coefficient (often an integer), default=0 [-]
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
    A : float
        Numerator coefficient, no default [varies]
    B : float
        Temperature exponent, no default [-]
    C : float, optional
        First denominator coefficient, default=0 [K]
    D : float, optional
        Second denominator coefficient, default=0 [K^2]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K;
                  property*K]

    Notes
    -----
    The derivative with respect to T is computed as follows.

    .. math::
        \frac{d Y}{dT} = \frac{A B T^{B}}{T \left(\frac{C}{T} + \frac{D}{T^{2}}
        + 1\right)} + \frac{A T^{B} \left(\frac{C}{T^{2}} + \frac{2 D}{T^{3}}
        \right)}{\left(\frac{C}{T} + \frac{D}{T^{2}} + 1\right)^{2}}


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
        if easy == 0.0:
            return easy
        try:
            return easy*T**B
        except:
            return 1e308
    elif order == 1:
        return (A*B*T**B/(T*(C/T + D/T**2 + 1))
                + A*T**B*(C/T**2 + 2*D/T**3)/(C/T + D/T**2 + 1)**2)
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
            x4 = x1*x3*trunc_log(B)
            r[0] = B**(-x2)
            r[1] = -x2*x3/B
            r[2] = -D*Ts[i]*x4/(C*C*x0)
            if x4 != 0.0:
                if x0 > 0:
                    r[3] = -x4*trunc_log(x0)
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
    A : float
        Constant coefficient [varies]
    B : float
        Temperature coefficient [K]
    C : float, optional
        Cubic temperature coefficient, default=0 [K^3]
    D : float, optional
        Power of 8 temperature coefficient, default=0 [K^8]
    E : float, optional
        Power of 9 temperature coefficient, default=0 [K^9]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -10, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == INTEGRAL_OVER_T_CALCULATION, unchanged from default]

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
    elif order == INTEGRAL_OVER_T_CALCULATION:
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
    A : float
        Multiplicative factor, [units]
    B : float
        Denominator power, [-]
    C : float
        Temperature denominator, [K]
    D : float
        Exponent for 1 - T/Tc usually, [-]
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
        problematic2 = problematic**D
        if abs(problematic2.imag) > 0.0: # This check should be removable - unless D is imaginary
            problematic2 = 0.0
        ans = A*B**(-(1. + problematic2))
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
    A : float
        Multiplier, [various]
    B : float
        Tau exponent constant term, [-]
    C : float, optional
        Tau exponent linear term, [-]
    D : float, optional
        Tau exponent quadratic term, [-]
    E : float, optional
        Tau exponent cubic term, [-]
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
        if tau <= 0.0:
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

def EQ106_AB(T, Tc, val, der):
    r'''Calculate the coefficients `A` and `B` of the DIPPR Equation #106 using
    the value of the function and its derivative at a specific point.

    .. math::
        A = val \left(\frac{1}{Tc} \left(- T + Tc\right)\right)^{- \frac{der}{val} \left(T - Tc\right)}

    .. math::
        B = \frac{der}{val} \left(T - Tc\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    val : float
        Property value [constant-specific]
    der : float
        First temperature derivative of property value [constant-specific/K]

    Returns
    -------
    A : float
        Parameter for the equation [constant-specific]
    B : float
        Parameter for the equation [-]

    Notes
    -----

    Examples
    --------
    >>> val = EQ106(300, 647.096, A=0.17766, B=2.567)
    >>> der = EQ106(300, 647.096, A=0.17766, B=2.567, order=1)
    >>> EQ106_AB(300, 647.096, val, der)
    (0.17766, 2.567)

    '''
    """# Derived with:
    from sympy import *
    T, Tc, A, B, val, der = symbols('T, Tc, A, B, val, der')
    Tr = T/Tc
    expr = A*(1 - Tr)**B

    Eq0 = Eq(expr, val)
    Eq1 = Eq(diff(expr, T), der)
    s = solve([Eq0, Eq1], [A, B])
    """
    x0 = T - Tc
    x1 = der*x0/val
    A, B = val*(-x0/Tc)**(-x1), x1
    return (A, B)

def EQ106_ABC(T, Tc, val, der, der2):
    r'''Calculate the coefficients `A`, `B`, and `C` of the DIPPR Equation #106
    using, the value of the function and its first and second derivative at a
    specific point.

    .. math::
        A = val \left(\frac{1}{Tc} \left(- T + Tc\right)\right)^{\frac{1}{val^{2}
        \left(\log{\left (\frac{1}{Tc} \left(- T + Tc\right) \right )} + 2\right)}
         \left(T \left(\log{\left (\frac{1}{Tc} \left(- T + Tc\right) \right )}
        + 1\right) \left(- T der^{2} + T der_{2} val + Tc der^{2} - Tc der_{2}
        val + der val\right) - T \left(- T der^{2} + T der_{2} val + Tc der^{2}
        - Tc der_{2} val + der val\right) - Tc \left(- T der^{2} + T der_{2} val
        + Tc der^{2} - Tc der_{2} val + der val\right) \log{\left (\frac{1}{Tc}
        \left(- T + Tc\right) \right )} - der val \left(T - Tc\right)
        \left(\log{\left (\frac{1}{Tc} \left(- T + Tc\right) \right )}
        + 2\right)\right)}

    .. math::
        B = \frac{1}{val^{2} \left(\log{\left (\frac{1}{Tc} \left(- T + Tc\right)
        \right )} + 2\right)} \left(- T \left(\log{\left (\frac{1}{Tc}
        \left(- T + Tc\right) \right )} + 1\right) \left(- T der^{2} + T der_{2}
        val + Tc der^{2} - Tc der_{2} val + der val\right) + Tc \left(- T der^{2}
        + T der_{2} val + Tc der^{2} - Tc der_{2} val + der val\right)
        \log{\left (\frac{1}{Tc} \left(- T + Tc\right) \right )} + der val
        \left(T - Tc\right) \left(\log{\left (\frac{1}{Tc} \left(- T + Tc\right) \right )} + 2\right)\right)

    .. math::
        C = \frac{Tc \left(- T der^{2} + T der_{2} val + Tc der^{2} - Tc der_{2}
        val + der val\right)}{val^{2} \left(\log{\left (\frac{1}{Tc}
        \left(- T + Tc\right) \right )} + 2\right)}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    val : float
        Property value [constant-specific]
    der : float
        First temperature derivative of property value [constant-specific/K]
    der2 : float
        Second temperature derivative of property value [constant-specific/K^2]

    Returns
    -------
    A : float
        Parameter for the equation [constant-specific]
    B : float
        Parameter for the equation [-]
    C : float
        Parameter for the equation [-]

    Notes
    -----

    Examples
    --------
    >>> val = EQ106(300, 647.096, A=0.17766, B=2.567, C=-0.01)
    >>> der = EQ106(300, 647.096, A=0.17766, B=2.567, C=-0.01, order=1)
    >>> der2 = EQ106(300, 647.096, A=0.17766, B=2.567, C=-0.01, order=2)
    >>> EQ106_ABC(300, 647.096, val, der, der2)
    (0.17766, 2.567, -0.01)

    '''
    """# Broken in recent versions of SymPy, SymPy 1.1 is good
    from sympy import *
    T, Tc, A, B, C, val, der, der2 = symbols('T, Tc, A, B, C, val, der, der2')
    Tr = T/Tc
    expr = A*(1 - Tr)**(B + C*Tr)

    Eq0 = Eq(expr, val)
    Eq1 = Eq(diff(expr, T), der)
    Eq2 = Eq(diff(expr, T, 2), der2)
    s = solve([Eq0, Eq1, Eq2], [A, B, C])
    """
    x0 = T - Tc
    x1 = -x0/Tc
    x2 = log(x1)
    x3 = x2 + 2
    x4 = 1/(val*val*x3)
    x5 = der*val
    x6 = der2*val
    x7 = der*der
    x8 = T*x6 - T*x7 - Tc*x6 + Tc*x7 + x5
    x9 = T*x8
    x10 = Tc*x8
    x11 = x0*x3*x5 + x10*x2 - x9*(x2 + 1)
    A, B, C = val*x1**(-x4*(x11 + x9)), x11*x4, x10*x4
    return (A, B, C)


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
    A : float, optional
        Constant property term, [J/(mol*K)]
    B : float, optional
        First hyperbolic term multiplier, [J/(mol*K)]
    C : float, optional
        First hyperbolic temperature denominator, [K]
    D : float, optional
        Second hyperbolic term multiplier, [J/(mol*K)]
    E : float, optional
        Second hyperbolic temperature denominator, [K]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -10, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == INTEGRAL_OVER_T_CALCULATION, unchanged from default]

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
    elif order == INTEGRAL_OVER_T_CALCULATION:
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
    A : float
        First coefficient, [-]
    B : float
        Second coefficient, [-]
    C : float
        Third coefficient, [-]
    D : float
        Fourth coefficient, [-]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -10, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == INTEGRAL_OVER_T_CALCULATION, unchanged from default]

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
        return A*A/t + 1.0*B + t*(-2.0*A*C + t*(-1.0*A*D + t*(-(1.0/3.0)*C*C + t*(-0.5*C*D - 0.2*D*D*t))))
        # return (A*A/t + B - 2.*A*C*t - A*D*t*t - (1.0/3.0)*C*C*t**3.
        #         - 0.5*C*D*t**4 - 0.2*D*D*t**5)
    elif order == 1:
        t = 1.-T/Tc
        return (A*A/(t*t) + 2.0*A*C + t*(2*A*D + t*(C*C + t*(2*C*D + D*D*t))))/Tc
    elif order == -1:
        x0 = D*D
        x1 = 2.0*D
        x2 = C*C
        x3 = C*D
        x4 = 6.0*x0
        x5 = A*D
        x6 = A*C
        T2 = T*T
        T3 = T2*T
        Tc2 = Tc*Tc
        Tc3 = Tc2*Tc
        return (-A*A*Tc*log(abs(T - Tc)) - D*T2*T3*(C + x1)/(10.0*Tc2*Tc2) + T3*T3*x0/(30.0*Tc2*Tc3) 
                + T2*T2*(x2 + 6.0*x3 + x4)/(12.0*Tc3) - T3*(2.0*x0 + x2 + 3.0*x3 + x5)/(3.0*Tc2) 
                + T2*(C*x1 + x0 + x2 + 2.0*x5 + 2.0*x6)/(2.0*Tc) 
                - T*(-30.0*B + 10.0*x2 + 15.0*x3 + x4 + 30.0*x5 + 60.0*x6)*(1.0/30))
    elif order == INTEGRAL_OVER_T_CALCULATION:
        x0 = A*A
        x1 = D*D
        x2 = 2.0*D
        x3 = C*C
        x4 = C*D
        x5 = 6.0*x1
        x6 = A*D
        x7 = A*C
        T2 = T*T
        T3 = T2*T
        Tc2 = Tc*Tc
        Tc3 = Tc2*Tc
        return (-D*T2*T2*(C + x2)/(8.0*Tc2*Tc2)
                 + T2*T3*x1/(25.0*Tc2*Tc3) + T3*(x3 + 6.0*x4 + x5)/(9.0*Tc3) 
                 - T2*(2.0*x1 + x3 + 3.0*x4 + x6)/(2.0*Tc2) 
                 + T*(C*x2 + x1 + x3 + 2.0*x6 + 2.0*x7)/Tc 
                 - x0*log(abs(T - Tc)) 
                 - (-30.0*B - 30.0*x0 + 10.0*x3 + 15.0*x4 + x5 + 30.0*x6 + 60.0*x7)*log(T)*(1.0/30.0))
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
    A : float
        First coefficient, [-]
    B : float
        Second coefficient, [K]
    C : float
        Third coefficient, [-]
    D : float
        Fourth coefficient, [1/K^2]
    E : float
        Fifth coefficient, [K^2]
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
    A : float
        First coefficient, [units]
    B : float
        Second coefficient, [units]
    C : float
        Third coefficient, [units]
    D : float
        Fourth coefficient, [units]
    E : float
        Fifth coefficient, [units]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -10, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == INTEGRAL_OVER_T_CALCULATION, unchanged from default]

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
    if T > Tc:
        T = Tc
    tau = 1.0-T/Tc
    cbrt_tau = cbrt(tau)
    if order == 0:
        return A + B*tau**0.35 + D*tau + C*cbrt_tau*cbrt_tau + E*tau*cbrt_tau
    elif order == 1:
        return (-0.35*B/((tau)**(0.65))
                - (2.0/3.0)*C/(cbrt_tau)
                - D - (4.0/3.0)*E*cbrt_tau)/Tc 
    elif order == -1:
        cbrt_tau2 = cbrt_tau*cbrt_tau
        cbrt_tau3 = cbrt_tau*cbrt_tau2
        return (A*T - (20.0/27)*B*Tc*(tau)**(1.35)
               + D*(-T*T/(2.0*Tc) + T)
               + cbrt_tau3*cbrt_tau2*(- 3.0/5.0*C*Tc
                - 3.0/7.0*E*Tc*cbrt_tau2))
    elif order == INTEGRAL_OVER_T_CALCULATION:
        # 3x increase in speed - cse via sympy
        x0 = log(T)
        x1 = 0.5*x0
        x2 = 1.0/Tc
        x3 = T*x2
        x4 = -x3 + 1.0
        x5 = 1.5*C
        x6 = cbrt(x4)
        x7 = 2*B
        x8 = x4**0.05
        x9 = log(-x6 + 1.0)
        x10 = 1.7320508075688772
        x11 = x10*atan(x10*((2/3.0)*x6 + 1.0/3.0))
        x12 = 2.23606797749979
        x13 = 0.5*x12
        x14 = x13 + 0.5
        x15 = B*x14
        x16 = sqrt(x13 + 2.5)
        x17 = 2.0*x8
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
        x28 = 3.1622776601683795*B/sqrt(x12 + 5.0)
        x29 = 2.0*x12
        x30 = sqrt(x29 + 10.0)
        x31 = 1.0/x30
        x32 = 1.0 - x12 
        x33 = 0.5*B*x22
        x34 = -x2*(T - Tc)
        x35 = 2.0*x34**0.1
        x36 = x35 + 2.0
        x37 = x34**0.05
        x38 = x30*x37
        x39 = 0.5*B*x16
        x40 = x37*sqrt(-x29 + 10.0)
        x41 = 0.25*x12
        x42 = B*(-x41 + 0.25)
        x43 = x12*x37
        x44 = x35 + x37 + 2.0
        x45 = B*(x41 + 0.25)
        x46 = -x43
        x47 = x35 - x37 + 2.0
        return A*x0 + 2.85714285714286*B*x4**0.35 - C*x1 + C*x11 + D*x0 - D*x3 - E*x1 - E*x11 + 0.75*E*x4**1.33333333333333 + 3.0*E*x6 + 1.5*E*x9 - x15*atan(x14*(x16 + x17)) + x15*atan(x14*(x16 + x18)) - x21*atan(x20*(x17 + x22)) + x21*atan(x20*(x18 + x22)) + x23*atan(x24*(x25 + x26)) - x23*atan(x24*(x25 + x27)) - x28*atan(x31*(x26 + x32)) + x28*atan(x31*(x27 + x32)) - x33*log(x36 - x38) + x33*log(x36 + x38) + x39*log(x36 - x40) - x39*log(x36 + x40) + x4**0.666666666666667*x5 - x42*log(x43 + x44) + x42*log(x46 + x47) + x45*log(x43 + x47) - x45*log(x44 + x46) + x5*x9 + x7*atan(x8) - x7*atanh(x8)
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
    A : float
        Constant property term, [J/(mol*K)]
    B : float
        First exponential term multiplier, [J/(mol*K)]
    C : float
        First exponential temperature denominator, [K]
    D : float
        Second exponential term multiplier, [J/(mol*K)]
    E : float
        Second exponential temperature denominator, [K]
    F : float
        Third exponential term multiplier, [J/(mol*K)]
    G : float
        Third exponential temperature denominator, [K]
    order : int, optional
        Order of the calculation. 0 for the calculation of the result itself;
        for 1, the first derivative of the property is returned, for
        -1, the indefinite integral of the property with respect to temperature
        is returned; and for -10, the indefinite integral of the property
        divided by temperature with respect to temperature is returned. No
        other integrals or derivatives are implemented, and an exception will
        be raised if any other order is given.

    Returns
    -------
    Y : float
        Property [constant-specific; if order == 1, property/K; if order == -1,
                  property*K; if order == INTEGRAL_OVER_T_CALCULATION, unchanged from default]

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
        T_inv = 1.0/T
        x0 = T_inv*T_inv
        x2 = exp(C*T_inv)
        x3 = exp(E*T_inv)
        x4 = exp(G*T_inv)
        x5 = x2 - 1.0
        x6 = x3 - 1.0
        x7 = x4 - 1.0
        return A + B*C*C*x0*x2/(x5*x5) + D*E*E*x0*x3/(x6*x6) + F*G*G*x0*x4/(x7*x7)
    elif order == 1:
        x0 = 1/T
        x1 = C*x0
        x2 = exp(x1)
        x3 = 1.0/(x2 - 1)
        x4 = x2*x3*x3
        x5 = E*x0
        x6 = exp(x5)
        x7 = 1.0/(x6 - 1)
        x8 = x6*x7*x7
        x9 = G*x0
        x10 = exp(x9)
        x11 = 1.0/(x10 - 1)
        x12 = x10*x11*x11
        x13 = C*C*C
        x14 = E*E*E
        x15 = G*G*G
        return (-2.0*B*C*C*x4 - B*x0*x13*x4 + 2.0*B*x0*x13*exp(2.0*x1)*x3*x3*x3 - 2.0*D*E*E*x8
                 - D*x0*x14*x8 + 2.0*D*x0*x14*exp(2.0*x5)*x7*x7*x7 - 2.0*F*G*G*x12 
                 - F*x0*x12*x15 + 2.0*F*x0*x15*exp(2.0*x9)*x11*x11*x11)*x0*x0*x0
    elif order == -1:
        T_inv = 1.0/T
        return (A*T + B*C*C/(C*exp(C*T_inv) - C) + D*E*E/(E*exp(E*T_inv) - E)
                + F*G*G/(G*exp(G*T_inv) - G))
    elif order == INTEGRAL_OVER_T_CALCULATION:
        x0 = 1.0/T
        x1 = exp(C*x0) - 1.0
        x2 = exp(E*x0) - 1.0
        x3 = exp(G*x0) - 1.0
        return A*log(T) + B*C*(x0 + x0/x1 - log(x1)/C) + D*E*(x0 + x0/x2 - log(x2)/E) + F*G*(x0 + x0/x3 - log(x3)/G)
    else:
        raise ValueError(order_not_found_msg)


dippr_eq_supported_orders = {
EQ100: (0, 1, -1, INTEGRAL_OVER_T_CALCULATION),
EQ101: (0, 1, 2, 3),
EQ102: (0, 1),
EQ104: (0, 1, -1, INTEGRAL_OVER_T_CALCULATION),
EQ105: (0, 1, 2, 3),
EQ106: (0, 1, 2, 3),
EQ107: (0, 1, -1, INTEGRAL_OVER_T_CALCULATION),
EQ114: (0, 1, -1, INTEGRAL_OVER_T_CALCULATION),
EQ115: (0, 1, 2, 3),
EQ116: (0, 1, -1, INTEGRAL_OVER_T_CALCULATION),
EQ127: (0, 1, -1, INTEGRAL_OVER_T_CALCULATION),
}
