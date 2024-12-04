"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021, 2022 Caleb Bell
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
.. autofunction:: chemicals.virial.BVirial_mixture
.. autofunction:: chemicals.virial.dBVirial_mixture_dzs
.. autofunction:: chemicals.virial.d2BVirial_mixture_dzizjs
.. autofunction:: chemicals.virial.d3BVirial_mixture_dzizjzks
.. autofunction:: chemicals.virial.CVirial_mixture_Orentlicher_Prausnitz
.. autofunction:: chemicals.virial.dCVirial_mixture_dT_Orentlicher_Prausnitz
.. autofunction:: chemicals.virial.d2CVirial_mixture_dT2_Orentlicher_Prausnitz
.. autofunction:: chemicals.virial.d3CVirial_mixture_dT3_Orentlicher_Prausnitz
.. autofunction:: chemicals.virial.dCVirial_mixture_Orentlicher_Prausnitz_dzs
.. autofunction:: chemicals.virial.d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs
.. autofunction:: chemicals.virial.d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks
.. autofunction:: chemicals.virial.d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs
.. autofunction:: chemicals.virial.dV_dzs_virial
.. autofunction:: chemicals.virial.d2V_dzizjs_virial


Second Virial Correlations
--------------------------
.. autofunction:: chemicals.virial.BVirial_Pitzer_Curl
.. autofunction:: chemicals.virial.BVirial_Abbott
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_extended

New implementations, returning the derivatives as well

.. autofunction:: chemicals.virial.BVirial_Pitzer_Curl_fast
.. autofunction:: chemicals.virial.BVirial_Abbott_fast
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_fast
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_extended_fast
.. autofunction:: chemicals.virial.BVirial_Oconnell_Prausnitz
.. autofunction:: chemicals.virial.BVirial_Xiang
.. autofunction:: chemicals.virial.BVirial_Meng
.. autofunction:: chemicals.virial.Meng_virial_a

Third Virial Correlations
-------------------------
.. autofunction:: chemicals.virial.CVirial_Orbey_Vera
.. autofunction:: chemicals.virial.CVirial_Liu_Xiang

Cross-Parameters
----------------
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_kijs
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_Tcijs
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_Pcijs
.. autofunction:: chemicals.virial.Tarakad_Danner_virial_CSP_omegaijs
.. autofunction:: chemicals.virial.Lee_Kesler_virial_CSP_Vcijs
.. autofunction:: chemicals.virial.Meng_Duan_2005_virial_CSP_kijs



Second Virial Correlations Dense Implementations
------------------------------------------------
.. autofunction:: chemicals.virial.BVirial_Xiang_vec
.. autofunction:: chemicals.virial.BVirial_Xiang_mat
.. autofunction:: chemicals.virial.BVirial_Pitzer_Curl_vec
.. autofunction:: chemicals.virial.BVirial_Pitzer_Curl_mat
.. autofunction:: chemicals.virial.BVirial_Abbott_vec
.. autofunction:: chemicals.virial.BVirial_Abbott_mat
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_vec
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_mat
.. autofunction:: chemicals.virial.BVirial_Meng_vec
.. autofunction:: chemicals.virial.BVirial_Meng_mat
.. autofunction:: chemicals.virial.BVirial_Oconnell_Prausnitz_vec
.. autofunction:: chemicals.virial.BVirial_Oconnell_Prausnitz_mat
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_extended_vec
.. autofunction:: chemicals.virial.BVirial_Tsonopoulos_extended_mat

Third Virial Correlations Dense Implementations
-----------------------------------------------
.. autofunction:: chemicals.virial.CVirial_Liu_Xiang_vec
.. autofunction:: chemicals.virial.CVirial_Orbey_Vera_vec
.. autofunction:: chemicals.virial.CVirial_Liu_Xiang_mat
.. autofunction:: chemicals.virial.CVirial_Orbey_Vera_mat
"""


__all__ = ['BVirial_Pitzer_Curl', 'BVirial_Pitzer_Curl_fast',
           'BVirial_Pitzer_Curl_vec', 'BVirial_Pitzer_Curl_mat',

           'BVirial_Abbott', 'BVirial_Abbott_fast',
           'BVirial_Abbott_vec', 'BVirial_Abbott_mat',

           'BVirial_Tsonopoulos', 'BVirial_Tsonopoulos_fast',
           'BVirial_Tsonopoulos_vec', 'BVirial_Tsonopoulos_mat',

           'BVirial_Tsonopoulos_extended',
           'BVirial_Tsonopoulos_extended_fast',
           'BVirial_Tsonopoulos_extended_vec',
           'BVirial_Tsonopoulos_extended_mat',

           'Meng_virial_a', 'BVirial_Meng',
           'BVirial_Meng_vec', 'BVirial_Meng_mat',

           'BVirial_Oconnell_Prausnitz','BVirial_Oconnell_Prausnitz_vec',
           'BVirial_Oconnell_Prausnitz_mat',

           'BVirial_Xiang', 'BVirial_Xiang_vec', 'BVirial_Xiang_mat',

           'BVirial_mixture', 'dBVirial_mixture_dzs',
           'd2BVirial_mixture_dzizjs', 'd3BVirial_mixture_dzizjzks',

           'dCVirial_mixture_Orentlicher_Prausnitz_dzs',
           'd2CVirial_mixture_Orentlicher_Prausnitz_dzizjs',
           'd3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks',

           'B_to_Z', 'B_from_Z', 'Z_from_virial_density_form',
           'Z_from_virial_pressure_form', 'CVirial_Orbey_Vera', 'CVirial_Liu_Xiang',
           'CVirial_Liu_Xiang_mat', 'CVirial_Liu_Xiang_vec',
           'CVirial_Orbey_Vera_vec', 'CVirial_Orbey_Vera_mat',
           'CVirial_mixture_Orentlicher_Prausnitz', 'dCVirial_mixture_dT_Orentlicher_Prausnitz',
           'd2CVirial_mixture_dT2_Orentlicher_Prausnitz',
           'd3CVirial_mixture_dT3_Orentlicher_Prausnitz',
           'd2CVirial_mixture_Orentlicher_Prausnitz_dTdzs',
           'Tarakad_Danner_virial_CSP_kijs', 'Tarakad_Danner_virial_CSP_Tcijs',
           'Tarakad_Danner_virial_CSP_Pcijs', 'Tarakad_Danner_virial_CSP_omegaijs',
           'Meng_Duan_2005_virial_CSP_kijs', 'Lee_Kesler_virial_CSP_Vcijs',
           'dV_dzs_virial', 'd2V_dzizjs_virial']


from fluids.constants import R, R_inv
from fluids.numerics import cbrt, exp, log, roots_cubic, sqrt
from fluids.numerics import numpy as np


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


def Z_from_virial_density_form(T, P, coeffs=None):
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
    coeffs : list[float], optional
        Sequence of virial coefficients [m³/mol^n], where n is the position:
        - First coefficient (n=1): Second virial coefficient B [m³/mol]
        - Second coefficient (n=2): Third virial coefficient C [m³/mol²]
        - Third coefficient (n=3): Fourth virial coefficient D [m³/mol³]
    
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
    >>> Z_from_virial_density_form(300, 122057.233762653, [1E-4, 1E-5, 1E-6, 1E-7])
    1.28434940526

    References
    ----------
    .. [1] Prausnitz, John M., Rudiger N. Lichtenthaler, and Edmundo Gomes de
       Azevedo. Molecular Thermodynamics of Fluid-Phase Equilibria. 3rd
       edition. Upper Saddle River, N.J: Prentice Hall, 1998.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering.
       Butterworth-Heinemann, 1985.
    '''
    if coeffs is None:
        return 1.0
    l = len(coeffs)
    if l == 1 or l==2 and coeffs[1] == 0.0:
        B = coeffs[0]
        determining_factor = 4.0*B*P + R*T
        if determining_factor < 0.0:
            Pmin = -0.25*R*T/(B)
#            raise ValueError("Too high of pressure for provided T and B") # numba: uncomment
            raise ValueError(f"Maximum allowed pressure for virial coefficient B of {B} m^3/mol is {Pmin} Pa") # numba: delete
        return 1.0/2. + 0.5*sqrt((determining_factor)/(R*T))
#        return ((R*T*(4*args[0]*P + R*T))**0.5 + R*T)/(2*P)
    if l == 2:
        B, C = coeffs[0], coeffs[1]
        sln = roots_cubic(-P/(R*T), 1.0, B, C)
        # print(sln)
        V = sln[0].real
        Z = P*V/(R*T)
        return Z
    size = l + 2
#    arr = np.ones(size, dtype=np.complex128) # numba: uncomment
    arr = [1.0]*size # numba: delete
    arr[-1] = -P/R/T
    for i in range(l):
        arr[-3-i] = coeffs[i]
    solns = np.roots(arr)
    for rho in solns:
        if abs(rho.imag) < 1e-12 and rho.real > 0.0:
            return float(P/(R*T*rho.real))
    raise ValueError("Could not find real root")


def Z_from_virial_pressure_form(P, coeffs=None):
    r'''Calculates the compressibility factor of a gas given its pressure, and
    pressure-form virial coefficients. Any number of coefficients is supported.

    .. math::
        Z = \frac{Pv}{RT} = 1 + B'P + C'P^2 + D'P^3 + E'P^4 \dots

    Parameters
    ----------
    P : float
        Pressure, [Pa]
    coeffs : list[float], optional
        List of pressure-form virial coefficients [1/Pa^n], where n is the index+1:
        coeffs[0]: Second virial coefficient B' [1/Pa]
        coeffs[1]: Third virial coefficient C' [1/Pa²]
        coeffs[2]: Fourth virial coefficient D' [1/Pa³]
        If None, returns ideal gas compressibility of 1.

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
    >>> Z_from_virial_pressure_form(102919.99946855308, [4.032286555169439e-09, 1.6197059494442215e-13, 6.483855042486911e-19])
    1.00283753944

    References
    ----------
    .. [1] Prausnitz, John M., Rudiger N. Lichtenthaler, and Edmundo Gomes de
       Azevedo. Molecular Thermodynamics of Fluid-Phase Equilibria. 3rd
       edition. Upper Saddle River, N.J: Prentice Hall, 1998.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering.
       Butterworth-Heinemann, 1985.
    '''
    if coeffs is None:
        return 1.0
    tot = 0.0
    fact = 1.0
    for i in range(len(coeffs)):
        tot += coeffs[i]*fact
        fact *= P
    return 1.0 + P*tot


### Second Virial Coefficients

def BVirial_mixture(zs, Bijs):
    r'''Calculate the `B` second virial coefficient from a matrix of
    virial cross-coefficients. The diagonal is virial coefficients of the
    pure components.

    .. math::
        B = \sum_i \sum_j y_i y_j B(T)

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Bijs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]

    Returns
    -------
    B : float
        Second virial coefficient in density form [m^3/mol]

    Notes
    -----

    Examples
    --------
    >>> Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    >>> zs = [.5, .3, .2]
    >>> BVirial_mixture(zs=zs, Bijs=Bijs)
    -3.19884e-05

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    B = 0.0
    N = len(Bijs)
    for i in range(N):
        B_tmp = 0.0
        row = Bijs[i]
        for j in range(N):
            B_tmp += zs[j]*row[j]
        B += zs[i]*(B_tmp)
    return B

def dBVirial_mixture_dzs(zs, Bijs, dB_dzs=None):
    r'''Calculate first mole fraction derivative of the `B` second virial
    coefficient from a matrix of virial cross-coefficients.

    .. math::
        \frac{\partial B}{\partial x_i} = \sum_j z_j(B_{i,j} + B_{j,i})

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Bijs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dzs : list[float], optional
        Array for first mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]

    Returns
    -------
    dB_dzs : list[float]
        First mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]

    Notes
    -----

    Examples
    --------
    >>> Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    >>> zs = [.5, .3, .2]
    >>> dBVirial_mixture_dzs(zs=zs, Bijs=Bijs)
    [-3.4089e-05, -7.2301e-05, -0.00012621]
    '''
    N = len(Bijs)
    if dB_dzs is None:
        dB_dzs = [0.0]*N
    else:
        for i in range(N):
            dB_dzs[i] = 0.0
    # Order this so that each row of B is processed sequentially.
    for k in range(N):
        zj = zs[k]
        Bks = Bijs[k]
        dB = 0.0
        for i in range(N):
            dB_dzs[i] += zj*Bks[i]
            # dB_dzs[k] +=zs[i]*Bks[i]
            dB +=zs[i]*Bks[i]
        dB_dzs[k] += dB

    # for k in range(N):
        # dB = 0.0
    #     Bks = Bijs[k]
    #     for i in range(N):
    #         dB += zs[i]*Bks[i]
    #     dB_dzs[k] += dB

    return dB_dzs

def d2BVirial_mixture_dzizjs(zs, Bijs, d2B_dzizjs=None):
    r'''Calculate second mole fraction derivative of the `B` second virial
    coefficient from a matrix of virial cross-coefficients.

    .. math::
        \frac{\partial^2 B}{\partial x_i \partial x_j} = B_{i,j} + B_{j,i}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Bijs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    d2B_dzizjs : list[list[float]], optional
        Array for First mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]


    Returns
    -------
    d2B_dzizjs : list[list[float]]
        First mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]

    Notes
    -----

    Examples
    --------
    >>> Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    >>> zs = [.5, .3, .2]
    >>> d2BVirial_mixture_dzizjs(zs=zs, Bijs=Bijs)
    [[-1.248e-05, -4.023e-05, -7.89e-05], [-4.023e-05, -8.782e-05, -0.0001292], [-7.89e-05, -0.0001292, -0.00024]]
    '''
    N = len(Bijs)
    if d2B_dzizjs is None:
        d2B_dzizjs = [[0.0]*N for _ in range(N)] # numba: delete
        # d2B_dzizjs = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Bi_row = Bijs[i]
        d2B_row = d2B_dzizjs[i]
        for j in range(N):
            d2B_row[j] = Bi_row[j] + Bijs[j][i]
    return d2B_dzizjs

def d3BVirial_mixture_dzizjzks(zs, Bijs, d3B_dzizjzks=None):
    r'''Calculate third mole fraction derivative of the `B` third virial
    coefficient from a matrix of virial cross-coefficients.

    .. math::
        \frac{\partial^3 B}{\partial x_i \partial x_j \partial x_k} = 0

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Bijs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    d3B_dzizjzks : list[list[list[float]]]
        Array for third mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]

    Returns
    -------
    d3B_dzizjzks : list[list[list[float]]]
        Third mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]

    Notes
    -----

    Examples
    --------
    >>> Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    >>> zs = [.5, .3, .2]
    >>> d3BVirial_mixture_dzizjzks(zs=zs, Bijs=Bijs)
    [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
    '''
    N = len(Bijs)
    if d3B_dzizjzks is None:
        d3B_dzizjzks = [[[0.0]*N for _ in range(N)] for _ in range(N)] # numba: delete
#        d3B_dzizjzks = np.zeros((N, N, N)) # numba: uncomment
    return d3B_dzizjzks


### B correlations

def BVirial_Oconnell_Prausnitz(T, Tc, Pc, omega):
    r'''Calculates the second virial coefficient using the model in [1]_.

    .. math::
        B_r=B^{(0)}+\omega B^{(1)}

    .. math::
        B^{(0)}=c0 + \frac{c1}{T_r} + \frac{c2}{T_r^2} + \frac{c3}{T_r^3}

    .. math::
        B^{(1)}=d0 + \frac{d1}{T_r^2} + \frac{d2}{T_r^3} + \frac{d3}{T_r^8}


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
    The coefficients are as follows:

    c0 = 0.1445
    c1 = -0.330
    c2 = -0.1385
    c3 = -0.0121

    d0 = 0.073
    d1 = 0.46
    d2 = -0.50
    d3 = -0.097
    d4 = -0.0073


    Examples
    --------
    >>> BVirial_Oconnell_Prausnitz(510., 425.2, 38E5, 0.193)
    (-0.000203193781, 1.036185972e-06, -6.53679132e-09, 6.59478287e-11)

    References
    ----------
    .. [1] O`Connell, J. P., and J. M. Prausnitz. "Empirical Correlation of
       Second Virial Coefficients for Vapor-Liquid Equilibrium Calculations."
       Industrial & Engineering Chemistry Process Design and Development 6,
       no. 2 (April 1, 1967): 245-50. https://doi.org/10.1021/i260022a016.
    '''
    c0 = 0.1445
    c1 = -0.330
    c2 = -0.1385
    c3 = -0.0121

    d0 = 0.073
    d1 = 0.46
    d2 = -0.50
    d3 = -0.097
    T_inv = 1.0/T
    T_inv2 = T_inv*T_inv
    x0 = Tc*T_inv
    T_inv3 = T_inv*T_inv2
    Tc2 = Tc*Tc
    Tc3 = Tc2*Tc
    x2 = Tc3*T_inv3
    x4 = Tc2*T_inv2
    x5 = R/Pc
    x6 = c2*x0
    x7 = c3*x4
    x8 = 2.0*d1
    x9 = d2*x0
    x10 = Tc3*Tc3*d3*T_inv3*T_inv3
    x11 = omega*x0
    x12 = Tc2*x5

    B = Tc*x5*(c0 + c1*x0 + c2*x4 + c3*x2 + omega*(d0 + d1*x4 + d2*x2 + Tc3*Tc3*Tc2*d3*T_inv3*T_inv3*T_inv2))
    dB = -x4*x5*(c1 + x11*(8.0*x10 + x8 + 3.0*x9) + 2.0*x6 + 3.0*x7)
    d2B = 2.0*T_inv3*x12*(c1 + 3.0*x11*(d1 + 12.0*x10 + 2.0*x9) + 3.0*x6 + 6.0*x7)
    d3B = -6.0*x12*(c1 + 2.0*x11*(60.0*x10 + x8 + 5.0*x9) + 4.0*x6 + 10.0*x7)*T_inv3*T_inv
    return (B, dB, d2B, d3B)

def BVirial_Oconnell_Prausnitz_vec(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the O'connell Prausnitz B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Oconnell_Prausnitz(T, Tcs[i], Pcs[i], omegas[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Oconnell_Prausnitz_mat(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Oconnell_Prausnitz B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        omega_row = omegas[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Oconnell_Prausnitz(T, Tc_row[j], Pc_row[j], omega_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

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

def BVirial_Pitzer_Curl_fast(T, Tc, Pc, omega):
    r'''Implementation of :obj:`BVirial_Pitzer_Curl` in the interface
    which calculates virial coefficients and their derivatives at the
    same time.

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
    >>> BVirial_Pitzer_Curl_fast(510., 425.2, 38E5, 0.193)
    (-0.000208453624, 1.065377516e-06, -5.7957101e-09, 4.513533043e-11)
    '''
    c0 = 0.1445
    c1 =  - 0.33
    c2 = - 0.1385
    c3 = - 0.0121

    d0 = 0.073
    d1 =  0.46
    d2 = - 0.5
    d3 = - 0.097
    d4 = - 0.0073

    T_inv = 1.0/T
    x0 = Tc*T_inv
    T_inv2 = T_inv*T_inv
    T_inv3 = T_inv2*T_inv
    T_inv4 = T_inv2*T_inv2
    Tc2 = Tc*Tc
    Tc3 = Tc*Tc2
    x2 = Tc3*T_inv3
    x4 = Tc2*T_inv2
    x5 = R/Pc
    x6 = 2.0*x0
    x7 = 3.0*x4
    x8 = Tc3*Tc3*Tc*d4*T_inv4*T_inv3
    x9 = 3.0*x0
    x10 = 6.0*x4
    x11 = Tc2*x5
    x12 = 4.0*x0
    x13 = 10.0*x4

    B = Tc*x5*(c0 + c1*x0 + c2*x4 + c3*x2 + omega*(d0 + d1*x0 + d2*x4 + d3*x2 + Tc3*Tc3*Tc2*d4*T_inv4*T_inv4))
    dB = -x4*x5*(c1 + c2*x6 + c3*x7 + omega*(d1 + d2*x6 + d3*x7 + 8.0*x8))
    d2B = 2.0*T_inv3*x11*(c1 + c2*x9 + c3*x10 + omega*(d1 + d2*x9 + d3*x10 + 36.0*x8))
    d3B = -6.0*x11*(c1 + c2*x12 + c3*x13 + omega*(d1 + d2*x12 + d3*x13 + 120.0*x8))*T_inv4
    return (B, dB, d2B, d3B)

def BVirial_Pitzer_Curl_vec(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the Pitzer-Curl B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Pitzer_Curl_fast(T, Tcs[i], Pcs[i], omegas[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Pitzer_Curl_mat(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Pitzer-Curl B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        omega_row = omegas[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Pitzer_Curl_fast(T, Tc_row[j], Pc_row[j], omega_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

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
    -0.000205701850095

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

def BVirial_Abbott_fast(T, Tc, Pc, omega):
    r'''Implementation of :obj:`BVirial_Abbott` in the interface
    which calculates virial coefficients and their derivatives at the
    same time.

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
    >>> BVirial_Abbott_fast(510., 425.2, 38E5, 0.193)
    (-0.0002057018500, 1.039249294e-06, -5.902233639e-09, 4.78222764e-11)
    '''
    c0 = 0.083
    c1 = -0.422
    d0 = 0.139
    d1 = - 0.172

    x0 = T/Tc
    T_inv = 1.0/T
    x1 = c1*x0**(-1.6)
    x2 = d1*x0**(-4.2)
    x3 = T_inv*R*Tc/Pc
    x4 = omega*x2

    B = x3*(c0 + omega*(d0 + x2) + x1)*T
    dB = -x3*(1.6*x1 + 4.2*x4)
    d2B = x3*(4.16*x1 + 21.84*x4)*T_inv
    d3B = -x3*(14.976*x1 + 135.408*x4)*T_inv*T_inv

    return (B, dB, d2B, d3B)

def BVirial_Abbott_vec(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the Abbott B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Abbott_fast(T, Tcs[i], Pcs[i], omegas[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Abbott_mat(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Abbott B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        omega_row = omegas[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Abbott_fast(T, Tc_row[j], Pc_row[j], omega_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

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
    -0.0002093529540

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

def BVirial_Tsonopoulos_fast(T, Tc, Pc, omega):
    r'''Implementation of :obj:`BVirial_Tsonopoulos` in the interface
    which calculates virial coefficients and their derivatives at the
    same time.

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
    >>> BVirial_Tsonopoulos_fast(510., 425.2, 38E5, 0.193)
    (-0.0002093529540, 9.95742355e-07, -5.54234465e-09, 4.57035160e-11)
    '''
    c0 = 0.1445
    c1 = -0.33
    c2 = -0.1385
    c3 = -0.0121
    c4 = -0.000607

    d0 = 0.0637
    d1 = 0.331
    d2 = - 0.423
    d3 = - 0.008

    T_inv = 1.0/T
    T_inv2 = T_inv*T_inv
    T_inv3 = T_inv2*T_inv
    T_inv4 = T_inv2*T_inv2
    Tc2 = Tc*Tc
    Tc4 = Tc2*Tc2
    x0 = Tc*T_inv
    x1 = Tc4*Tc4*T_inv4*T_inv4
    x3 = Tc2*Tc*T_inv3
    x5 = Tc2*T_inv2
    x6 = R/Pc
    x7 = c2*x0
    x8 = Tc4*Tc2*Tc*c4*T_inv3*T_inv4
    x9 = c3*x5
    x10 = 2.0*d1
    x11 = d2*x0
    x12 = Tc4*Tc2*d3*T_inv3*T_inv3
    x13 = omega*x0
    x14 = Tc2*x6


    B = Tc*x6*(c0 + c1*x0 + c2*x5 + c3*x3 + c4*x1 + omega*(d0 + d1*x5 + d2*x3 + d3*x1))
    dB = -x5*x6*(c1 + x13*(x10 + 3.0*x11 + 8.0*x12) + 2.0*x7 + 8.0*x8 + 3.0*x9)
    d2B = 2.0*x14*T_inv3*(c1 + 3.0*x13*(d1 + 2.0*x11 + 12.0*x12) + 3.0*x7 + 36.0*x8 + 6.0*x9)
    d3B = -6.0*x14*(c1 + 2.0*x13*(x10 + 5.0*x11 + 60.0*x12) + 4.0*x7 + 120.0*x8 + 10.0*x9)*T_inv4

    return (B, dB, d2B, d3B)

def BVirial_Tsonopoulos_vec(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the Tsonopoulos B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Tsonopoulos_fast(T, Tcs[i], Pcs[i], omegas[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Tsonopoulos_mat(T, Tcs, Pcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Tsonopoulos B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        omega_row = omegas[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Tsonopoulos_fast(T, Tc_row[j], Pc_row[j], omega_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

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
    -9.679718337596e-05

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
        if species_type in ('simple', 'normal'):
            a, b = 0, 0
        elif species_type == 'methyl alcohol':
            a, b = 0.0878, 0.0525
        elif species_type == 'water':
            a, b = -0.0109, 0
        elif dipole != 0 and Tc != 0 and Pc != 0:
            dipole_r = 1E5*dipole**2*(Pc/101325.0)/Tc**2

            if species_type in ('ketone', 'aldehyde', 'alkyl nitrile',  'ether', 'carboxylic acid', 'ester'):
                a, b = -2.14E-4*dipole_r-4.308E-21*dipole_r**8, 0.0
            elif species_type in ('alkyl halide', 'mercaptan', 'sulfide', 'disulfide'):
                a, b = -2.188E-4*dipole_r**4-7.831E-21*dipole_r**8, 0.0

            elif species_type == 'alkanol':
                a, b = 0.0878, 0.00908+0.0006957*dipole_r
    Br = B0 + omega*B1 + a*B2 + b*B3
    return Br*R*Tc/Pc

def BVirial_Tsonopoulos_extended_fast(T, Tc, Pc, omega, a=0.0, b=0.0):
    r'''Implementation of :obj:`BVirial_Tsonopoulos_extended` in the interface
    which calculates virial coefficients and their derivatives at the
    same time.

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
        Fit parameter [-]
    b : float, optional
        Fit parameter [-]

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
    >>> BVirial_Tsonopoulos_extended_fast(510., 425.2, 38E5, 0.193)
    (-0.0002093529540, 9.9574235e-07, -5.54234465e-09, 4.5703516e-11)
    '''
    c0 = 0.1445
    c1 = -0.33
    c2 = -0.1385
    c3 = -0.0121
    c4 = -0.000607

    d0 = 0.0637
    d1 = 0.331
    d2 = -0.423
    d3 = -0.008

    T_inv = 1.0/T
    T_inv2 = T_inv*T_inv
    T_inv3 = T_inv*T_inv2
    T_inv4 = T_inv2*T_inv2
    Tc2 = Tc*Tc
    Tc3 = Tc*Tc2
    Tc4 = Tc*Tc3
    x0 = Tc*T_inv
    x1 = Tc4*Tc4*T_inv4*T_inv4
    x2 = Tc3*Tc3*T_inv3*T_inv3
    x4 = Tc3*T_inv3
    x6 = Tc2*T_inv2
    x7 = R/Pc
    x8 = c2*x0
    x9 = Tc4*Tc3*T_inv3*T_inv4
    x10 = 8.0*x9
    x11 = Tc3*Tc2*a*T_inv3*T_inv2
    x12 = c3*x6
    x13 = 2.0*d1
    x14 = d2*x0
    x15 = d3*x2
    x16 = omega*x0
    x17 = 36.0*x9
    x18 = Tc2*x7
    x19 = 120.0*x9
    B = Tc*x7*(a*x2 - b*x1 + c0 + c1*x0 + c2*x6 + c3*x4 + c4*x1 + omega*(d0 + d1*x6 + d2*x4 + d3*x1))
    dB = -x6*x7*(-b*x10 + c1 + c4*x10 + 6.0*x11 + 3.0*x12 + x16*(x13 + 3.0*x14 + 8.0*x15) + 2.0*x8)
    d2B = 2.0*x18*T_inv3*(-b*x17 + c1 + c4*x17 + 21.0*x11 + 6.0*x12 + 3.0*x16*(d1 + 2.0*x14 + 12.0*x15) + 3.0*x8)
    d3B = -6.0*x18*(-b*x19 + c1 + c4*x19 + 56.0*x11 + 10.0*x12 + 2.0*x16*(x13 + 5.0*x14 + 60.0*x15) + 4.0*x8)*T_inv4
    return (B, dB, d2B, d3B)


def BVirial_Tsonopoulos_extended_vec(T, Tcs, Pcs, omegas, ais, bs, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the Tsonopoulos (extended) B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    omegas : list[float]
        Acentric factor for fluids, [-]
    ais : list[float]
        Fit parameters, [-]
    bs : list[float]
        Fit parameters, [-]
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Tsonopoulos_extended_fast(T, Tcs[i], Pcs[i], omegas[i], ais[i], bs[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Tsonopoulos_extended_mat(T, Tcs, Pcs, omegas, ais, bs, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Tsonopoulos (extended) B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    ais : list[list[float]]
        Fit parameters, [-]
    bs : list[list[float]]
        Fit parameters, [-]
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        omega_row = omegas[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]
        a_row = ais[i]
        b_row = bs[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Tsonopoulos_extended_fast(T, Tc_row[j], Pc_row[j], omega_row[j], a_row[j], b_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s


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

def BVirial_Xiang_vec(T, Tcs, Pcs, Vcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the Xiang B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    Vcs : list[float]
        Critical volume of the fluids [m^3/mol]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Xiang(T, Tcs[i], Pcs[i], Vcs[i], omegas[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Xiang_mat(T, Tcs, Pcs, Vcs, omegas, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Xiang B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    Vcs : list[list[float]]
        Critical volume of the fluids [m^3/mol]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        Vc_row = Vcs[i]
        omega_row = omegas[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Xiang(T, Tc_row[j], Pc_row[j], Vc_row[j], omega_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Meng(T, Tc, Pc, Vc, omega, a=0.0):
    r'''Calculates the second virial coefficient using the model in [1]_.

    .. math::
        B = \frac{RT_c}{P_c}\left(f_0 + \omega f_1 + f_2\right)

    .. math::
        f_0 = c_0 + c_1/T_r + c_2/T_r^2 + c_3/T_r^3 + c_4/T_r^8

    .. math::
        f_1 = d_0 + d_1/T_r + d_2/T_r^2 + d_3/T_r^3 + d_4/T_r^8

    .. math::
        f_2 = \frac{a}{T_r^6}

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
    a : float
        Polar parameter that can be estimated by :obj:`chemicals.virial.Meng_virial_a`

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
    >>> BVirial_Meng(388.26, 647.1, 22050000.0, 5.543076e-05, 0.344)
    (-0.00032436028, 2.47004e-06, -3.132e-08, 5.8e-10)

    References
    ----------
    .. [1] Meng, Long, Yuan-Yuan Duan, and Lei Li. "Correlations for Second and
       Third Virial Coefficients of Pure Fluids." Fluid Phase Equilibria 226
       (December 10, 2004): 109-20. https://doi.org/10.1016/j.fluid.2004.09.023.
    '''
    c0, c1, c2, c3, c4 =  0.13356, -0.30252, -0.15668, -0.00724, -0.00022
    d0, d1, d2, d3, d4 = 0.17404, -0.15581, 0.38183, -0.44044, -0.00541
    T_inv = 1.0/T
    Tc_T = Tc*T_inv
    Tc_T2 = Tc_T*Tc_T
    Tc_T4 = Tc_T2*Tc_T2
    x1 = Tc_T4*Tc_T4
    T_inv3 = T_inv*T_inv*T_inv
    Tc2 = Tc*Tc
    x3 = Tc*Tc2*T_inv3
    x5 = Tc2*T_inv*T_inv
    x6 = R/Pc
    x7 = 2.0*Tc_T
    x8 = Tc_T4*Tc_T2*Tc_T
    x9 = 8.0*x8
    x10 = a*Tc_T4*Tc_T
    x11 = 3.0*x5
    x12 = 3.0*Tc_T
    x13 = 36.0*x8
    x14 = 6.0*x5
    x15 = Tc2*x6
    x16 = 4.0*Tc_T
    x17 = 120.0*x8
    x18 = 10.0*x5
    B = Tc*x6*(c0 + c1*Tc_T + c2*x5 + c3*x3 + c4*x1 + omega*(d0 + d1*Tc_T + d2*x5 + d3*x3 + d4*x1) + a*Tc_T4*Tc_T2)
    dB = -x5*x6*(c1 + c2*x7 + c3*x11 + c4*x9 + omega*(d1 + d2*x7 + d3*x11 + d4*x9) + 6.0*x10)
    dB2 = 2.0*x15*T_inv3*(c1 + c2*x12 + c3*x14 + c4*x13 + omega*(d1 + d2*x12 + d3*x14 + d4*x13) + 21.0*x10)
    dB3 = -6.0*x15*(c1 + c2*x16 + c3*x18 + c4*x17 + omega*(d1 + d2*x16 + d3*x18 + d4*x17) + 56.0*x10)*T_inv3*T_inv
    return B, dB, dB2, dB3

def BVirial_Meng_vec(T, Tcs, Pcs, Vcs, omegas, ais, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a vectorized calculation of the Meng B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    Vcs : list[float]
        Critical volume of the fluids [m^3/mol]
    omegas : list[float]
        Acentric factor for fluids, [-]
    ais : list[float]
        Polar parameters that can be estimated by :obj:`chemicals.virial.Meng_virial_a`
    Bs : list[float], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[float]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[float]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[float]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[float]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [0.0]*N
    if dB_dTs is None:
        dB_dTs = [0.0]*N
    if d2B_dT2s is None:
        d2B_dT2s = [0.0]*N
    if d3B_dT3s is None:
        d3B_dT3s = [0.0]*N
    for i in range(N):
        B, dB, d2B, d3B = BVirial_Meng(T, Tcs[i], Pcs[i], Vcs[i], omegas[i], ais[i])
        Bs[i] = B
        dB_dTs[i] = dB
        d2B_dT2s[i] = d2B
        d3B_dT3s[i] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def BVirial_Meng_mat(T, Tcs, Pcs, Vcs, omegas, ais, Bs=None, dB_dTs=None,
                      d2B_dT2s=None, d3B_dT3s=None):
    r'''Perform a matrix calculation of the Meng B virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    Vcs : list[list[float]]
        Critical volume of the fluids [m^3/mol]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    ais : list[float]
        Polar parameters that can be estimated as the average of the pure
        component values predicted by :obj:`chemicals.virial.Meng_virial_a`
    Bs : list[list[float]], optional
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]], optional
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]], optional
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]], optional
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Returns
    -------
    Bs : list[list[float]]
        Second virial coefficient in density form [m^3/mol]
    dB_dTs : list[list[float]]
        First temperature derivative of second virial coefficient in density
        form [m^3/mol/K]
    d2B_dT2s : list[list[float]]
        Second temperature derivative of second virial coefficient in density
        form [m^3/mol/K^2]
    d3B_dT3s : list[list[float]]
        Third temperature derivative of second virial coefficient in density
        form [m^3/mol/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Bs is None:
        Bs = [[0.0]*N for _ in range(N)] # numba: delete
#        Bs = np.zeros((N, N)) # numba: uncomment
    if dB_dTs is None:
        dB_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dB_dTs = np.zeros((N, N)) # numba: uncomment
    if d2B_dT2s is None:
        d2B_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2B_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3B_dT3s is None:
        d3B_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3B_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        Vc_row = Vcs[i]
        omega_row = omegas[i]
        a_row = ais[i]

        B_row = Bs[i]
        dB_row = dB_dTs[i]
        d2B_row = d2B_dT2s[i]
        d3B_row = d3B_dT3s[i]

        for j in range(N):
            B, dB, d2B, d3B = BVirial_Meng(T, Tc_row[j], Pc_row[j], Vc_row[j], omega_row[j], a_row[j])
            B_row[j] = B
            dB_row[j] = dB
            d2B_row[j] = d2B
            d3B_row[j] = d3B
    return Bs, dB_dTs, d2B_dT2s, d3B_dT3s

def Meng_virial_a(Tc, Pc, dipole=0.0, haloalkane=False):
    r'''Calculate the `a` parameter which is used in the Meng
    `B` second virial coefficient for polar components. There are two
    correlations implemented - one for haloalkanes, and another for other
    polar molecules. If the dipole moment is not provided, a value of 0.0
    will be returned.

    If the compound is a haloalkane

    .. math::
        a = -1.1524\times 10^{-6}{\mu}_r^2 + 7.2238\times 10^{-11}{\mu}_r^4
        - 1.8701\times 10^{-15}{\mu}_r^6

    Otherwise

    .. math::
        a = -3.0309\times 10^{-6}{\mu}_r^2 + 9.503\times 10^{-11}{\mu}_r^4
        - 1.2469\times 10^{-15}{\mu}_r^6

    Parameters
    ----------
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    dipole : float
        Dipole moment, [debye]
    haloalkane : bool
        Whether or not the compound is a haloalkane, [-]

    Returns
    -------
    a : float
        Coefficient [-]

    Notes
    -----

    Examples
    --------
    Ethanol

    >>> Meng_virial_a(514.0, 6137000.0, 1.44, haloalkane=False)
    -0.00637841

    R-41 Fluoromethane

    >>> Meng_virial_a(317.4, 5870000.0, 1.85, haloalkane=True)
    -0.04493829

    References
    ----------
    .. [1] Meng, Long, Yuan-Yuan Duan, and Lei Li. "Correlations for Second and
       Third Virial Coefficients of Pure Fluids." Fluid Phase Equilibria 226
       (December 10, 2004): 109-20. https://doi.org/10.1016/j.fluid.2004.09.023.
    '''
    # Perfect validated with graph
    mur = dipole*dipole*Pc/(1.01325*Tc*Tc)
    if haloalkane:
        a = -1.1524e-6*mur**2 + 7.2238e-11*mur**4 - 1.8701E-15*mur**6
    else:
        a = -3.0309E-6*mur**2 + 9.503E-11*mur**4 - 1.2469E-15*mur**6
    return a

def Kronecker_delta(i, j): return 1 if i == j else 0.0


def dCVirial_mixture_Orentlicher_Prausnitz_dzs(zs, Cijs, dCs=None):
    r'''Calculate the first mole fraction derivatives of the `C` third virial
    coefficient from a matrix of
    virial cross-coefficients.

    .. math::
        \frac{\partial C}{\partial z_m} =
        \sum_{\substack{0 \leq i \leq nc\\0 \leq j \leq nc\\0 \leq k \leq nc}}
        \sqrt[3]{{Cs}_{i,j} {Cs}_{i,k} {Cs}_{j,k}} \left(\delta_{i m} {zs}_{j}
        {zs}_{k} + \delta_{j m} {zs}_{i} {zs}_{k} + \delta_{k m} {zs}_{i}
        {zs}_{j}\right)

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    dCs : list[float], optional
        First derivatives of C with respect to mole fraction, [m^6/mol^2]

    Returns
    -------
    dC_dzs : list[float]
        First derivatives of C with respect to mole fraction, [m^6/mol^2]

    Notes
    -----
    This equation can be derived with SymPy, as follows

    >>> from sympy import * # doctest: +SKIP
    >>> i, j, k, m, n, o = symbols("i, j, k, m, n, o", cls=Idx) # doctest: +SKIP
    >>> zs = IndexedBase('zs') # doctest: +SKIP
    >>> Cs = IndexedBase('Cs') # doctest: +SKIP
    >>> nc = symbols('nc') # doctest: +SKIP
    >>> C_expr = Sum(zs[i]*zs[j]*zs[k]*cbrt(Cs[i,j]*Cs[i,k]*Cs[j,k]),[i,0,nc],[j,0,nc],[k,0,nc]) # doctest: +SKIP
    >>> diff(C_expr, zs[m]) # doctest: +SKIP
    Sum((Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*KroneckerDelta(i, m)*zs[j]*zs[k] + (Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*KroneckerDelta(j, m)*zs[i]*zs[k] + (Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*KroneckerDelta(k, m)*zs[i]*zs[j], (i, 0, nc), (j, 0, nc), (k, 0, nc))

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> zs = [.5, .3, .2]
    >>> dCVirial_mixture_Orentlicher_Prausnitz_dzs(zs, Cijs)
    [5.44450470e-09, 6.54968776e-09, 7.74986672e-09]
    '''
    N = len(zs)
    Cij_cbrts = [[0.0]*N for _ in range(N)] # numba: delete
#    Cij_cbrts = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Cij_cbrt_row = Cij_cbrts[i]
        Cij_row = Cijs[i]
        for j in range(i):
            if Cij_row[j] > 0.0:
                Cij_cbrt_row[j] = Cij_cbrts[j][i] = cbrt(Cij_row[j])

        if Cij_row[i] > 0.0:
            Cij_cbrt_row[i] = cbrt(Cij_row[i])

    cC = Cij_cbrts

    if dCs is None:
        dCs = [0.0]*N


    for m in range(N):
        dC = 0.0
        for i in range(N):
            fact = cC[i][m]*zs[i]
            tmp = 0.0
            for j in range(N):
                tmp += cC[i][j]*cC[m][j]*zs[j]
            dC += fact*tmp

        # symmetry double this term
        dC += dC
        dCs[m] += dC
    for j in range(N):
        for k in range(N):
            x0 = zs[j]*zs[k]*cC[j][k]
            for m in range(N):
                dCs[m] += x0*cC[m][j]*cC[m][k]
    return dCs

def d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs(zs, Cijs, d2Cs=None):
    r'''Calculate the second mole fraction derivatives of the `C` third virial
    coefficient from a matrix of
    virial cross-coefficients.

    .. math::
        \frac{\partial^2 C}{\partial z_m \partial z_n} =
        \sum_{\substack{0 \leq i \leq nc\\0 \leq j \leq nc\\0 \leq k \leq nc}}
        \sqrt[3]{{Cs}_{i,j} {Cs}_{i,k} {Cs}_{j,k}} \left(\delta_{i m}
        \delta_{j n} {zs}_{k} + \delta_{i m} \delta_{k n} {zs}_{j}
        + \delta_{i n} \delta_{j m} {zs}_{k} + \delta_{i n} \delta_{k m}
        {zs}_{j} + \delta_{j m} \delta_{k n} {zs}_{i}
        + \delta_{j n} \delta_{k m} {zs}_{i}\right)

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    d2Cs : list[list[float]], optional
        Second derivatives of C with respect to mole fraction, [m^6/mol^2]

    Returns
    -------
    d2Cs : list[list[float]]
        Second derivatives of C with respect to mole fraction, [m^6/mol^2]

    Notes
    -----
    This equation can be derived with SymPy, as follows

    >>> from sympy import * # doctest: +SKIP
    >>> i, j, k, m, n, o = symbols("i, j, k, m, n, o", cls=Idx) # doctest: +SKIP
    >>> zs = IndexedBase('zs') # doctest: +SKIP
    >>> Cs = IndexedBase('Cs') # doctest: +SKIP
    >>> nc = symbols('nc') # doctest: +SKIP
    >>> C_expr = Sum(zs[i]*zs[j]*zs[k]*cbrt(Cs[i,j]*Cs[i,k]*Cs[j,k]),[i,0,nc],[j,0,nc],[k,0,nc]) # doctest: +SKIP
    >>> diff(C_expr, zs[m], zs[n]) # doctest: +SKIP
    Sum((Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*(KroneckerDelta(i, m)*KroneckerDelta(j, n)*zs[k] + KroneckerDelta(i, m)*KroneckerDelta(k, n)*zs[j] + KroneckerDelta(i, n)*KroneckerDelta(j, m)*zs[k] + KroneckerDelta(i, n)*KroneckerDelta(k, m)*zs[j] + KroneckerDelta(j, m)*KroneckerDelta(k, n)*zs[i] + KroneckerDelta(j, n)*KroneckerDelta(k, m)*zs[i]), (i, 0, nc), (j, 0, nc), (k, 0, nc))

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> zs = [.5, .3, .2]
    >>> d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs(zs, Cijs)
    [[9.6827886655e-09, 1.1449146725e-08, 1.3064355337e-08], [1.1449146725e-08, 1.38557674294e-08, 1.60903596751e-08], [1.3064355337e-08, 1.60903596751e-08, 2.0702239403e-08]]
    '''
    N = len(zs)
    Cij_cbrts = [[0.0]*N for _ in range(N)] # numba: delete
#    Cij_cbrts = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Cij_cbrt_row = Cij_cbrts[i]
        Cij_row = Cijs[i]
        for j in range(i):
            if Cij_row[j] > 0.0:
                Cij_cbrt_row[j] = Cij_cbrts[j][i] = Cij_row[j]**(1.0/3)

        if Cij_row[i] > 0.0:
            Cij_cbrt_row[i] = Cij_row[i]**(1.0/3.0)

    cC = Cij_cbrts

    if d2Cs is None:
        d2Cs = [[0.0]*N for _ in range(N)] # numba: delete
        # d2Cs = np.zeros((N, N)) # numba: uncomment

    for i in range(N):
        cCi = cC[i]
        zi = zs[i]
        d2Cis = d2Cs[i]
        for j in range(N):
            cCij = cC[i][j]
            cCj = cC[j]
            zj = zs[j]
            d2Cjs = d2Cs[j]
            for k in range(N):
                cCv = cCij*cCi[k]*cCj[k]

                d2Cis[j] += cCv*zs[k]
                d2Cis[k] += cCv*zj
                d2Cjs[i] += cCv*zs[k]
                d2Cjs[k] += cCv*zi
                d2Cs[k][i] += cCv*zj
                d2Cs[k][j] += cCv*zi
    return d2Cs


def d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks(zs, Cijs, d3Cs=None):
    r'''Calculate the third mole fraction derivatives of the `C` third virial
    coefficient from a matrix of
    virial cross-coefficients.

    .. math::
        \frac{\partial^3 C}{\partial z_m \partial z_n \partial z_o} =
        \sum_{\substack{0 \leq i \leq nc\\0 \leq j \leq nc\\0 \leq k \leq nc}}
        \sqrt[3]{{Cs}_{i,j} {Cs}_{i,k} {Cs}_{j,k}} \left(\delta_{i m} \delta_{j n}
        \delta_{k o} + \delta_{i m} \delta_{j o} \delta_{k n} + \delta_{i n}
        \delta_{j m} \delta_{k o} + \delta_{i n} \delta_{j o} \delta_{k m}
        + \delta_{i o} \delta_{j m} \delta_{k n} + \delta_{i o} \delta_{j n}
        \delta_{k m}\right)

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    d3Cs : list[list[list[float]]], optional
        Third derivatives of C with respect to mole fraction, [m^6/mol^2]

    Returns
    -------
    d3Cs : list[list[list[float]]]
        Third derivatives of C with respect to mole fraction, [m^6/mol^2]

    Notes
    -----
    This equation can be derived with SymPy, as follows

    >>> from sympy import * # doctest: +SKIP
    >>> i, j, k, m, n, o = symbols("i, j, k, m, n, o", cls=Idx) # doctest: +SKIP
    >>> zs = IndexedBase('zs') # doctest: +SKIP
    >>> Cs = IndexedBase('Cs') # doctest: +SKIP
    >>> nc = symbols('nc') # doctest: +SKIP
    >>> C_expr = Sum(zs[i]*zs[j]*zs[k]*cbrt(Cs[i,j]*Cs[i,k]*Cs[j,k]),[i,0,nc],[j,0,nc],[k,0,nc]) # doctest: +SKIP
    >>> diff(C_expr, zs[m], zs[n], zs[o]) # doctest: +SKIP
    Sum((Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*(KroneckerDelta(i, m)*KroneckerDelta(j, n)*KroneckerDelta(k, o) + KroneckerDelta(i, m)*KroneckerDelta(j, o)*KroneckerDelta(k, n) + KroneckerDelta(i, n)*KroneckerDelta(j, m)*KroneckerDelta(k, o) + KroneckerDelta(i, n)*KroneckerDelta(j, o)*KroneckerDelta(k, m) + KroneckerDelta(i, o)*KroneckerDelta(j, m)*KroneckerDelta(k, n) + KroneckerDelta(i, o)*KroneckerDelta(j, n)*KroneckerDelta(k, m)), (i, 0, nc), (j, 0, nc), (k, 0, nc))

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> zs = [.5, .3, .2]
    >>> d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks(zs, Cijs)
    [[[8.760000000e-09, 1.0187346981e-08, 1.12329228549e-08], [1.01873469818e-08, 1.21223973593e-08, 1.35937701316e-08], [1.12329228549e-08, 1.35937701316e-08, 1.68488143533e-08]], [[1.01873469818e-08, 1.21223973593e-08, 1.35937701316e-08], [1.2122397359e-08, 1.47600000000e-08, 1.68328437491e-08], [1.35937701316e-08, 1.68328437491e-08, 2.12181074230e-08]], [[1.12329228549e-08, 1.35937701316e-08, 1.68488143533e-08], [1.35937701316e-08, 1.68328437491e-08, 2.12181074230e-08], [1.68488143533e-08, 2.12181074230e-08, 2.9562000000e-08]]]
    '''
    N = len(zs)
    Cij_cbrts = [[0.0]*N for _ in range(N)] # numba: delete
#    Cij_cbrts = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Cij_cbrt_row = Cij_cbrts[i]
        Cij_row = Cijs[i]
        for j in range(i):
            if Cij_row[j] > 0.0:
                Cij_cbrt_row[j] = Cij_cbrts[j][i] = cbrt(Cij_row[j])

        if Cij_row[i] > 0.0:
            Cij_cbrt_row[i] = cbrt(Cij_row[i])

    cC = Cij_cbrts

    if d3Cs is None:
        d3Cs = [[[0.0]*N for _ in range(N)] for _ in range(N)]# numba: delete
        # d3Cs = np.zeros((N, N, N)) # numba: uncomment
    for i in range(N):
        d3Cis = d3Cs[i]
        cCis = cC[i]
        for j in range(N):
            d3Cijs = d3Cis[j]
            d3Cjs = d3Cs[j]
            cCij = cC[i][j]
            cCjs = cC[j]
            for k in range(N):
                cCv = cCij*cCis[k]*cCjs[k]
                d3Cijs[k]+= cCv
                d3Cis[k][j]+= cCv
                d3Cjs[i][k]+= cCv
                d3Cjs[k][i]+= cCv
                d3Cs[k][i][j]+= cCv
                d3Cs[k][j][i]+= cCv
    return d3Cs

def d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs(zs, Cijs, dCij_dTs,
                                                  d2C_dTdzs=None):
    r'''Calculate the first mole fraction derivatives of the `C` third virial
    coefficient from a matrix of
    virial cross-coefficients.

    .. math::
        \frac{\partial^2 C}{\partial T \partial z_m}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    dCij_dTs : list[list[float]]
        First temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K]
    d2C_dTdzs : list[float], optional
        Array for second derivatives of C with respect to mole fraction and
        temperature, [m^6/mol^2/K]

    Returns
    -------
    d2C_dTdzs : list[float]
        Second derivatives of C with respect to mole fraction and temperature,
        [m^6/mol^2/K]

    Notes
    -----
    This equation can be derived with SymPy, as follows

    >>> from sympy import * # doctest: +SKIP
    >>> from sympy import * # doctest: +SKIP
    >>> i, j, k, m, n, o, T = symbols("i, j, k, m, n, o, T", cls=Idx) # doctest: +SKIP
    >>> zs = IndexedBase('zs') # doctest: +SKIP
    >>> Cs = IndexedBase('Cs') # doctest: +SKIP
    >>> dC_dTs = IndexedBase('dC_dTs') # doctest: +SKIP
    >>> nc = symbols('nc') # doctest: +SKIP
    >>> C_expr = Sum(zs[i]*zs[j]*zs[k]/3*cbrt(Cs[i,j]*Cs[i,k]*Cs[j,k])/(Cs[i,j]*Cs[i,k]*Cs[j,k])*(Cs[i,j]*Cs[i,k]*dC_dTs[j,k] + Cs[i,j]*dC_dTs[i,k]*Cs[j,k] + dC_dTs[i,j]*Cs[i,k]*Cs[j,k]),[i,0,nc],[j,0,nc],[k,0,nc]) # doctest: +SKIP
    >>> diff(C_expr, zs[m]) # doctest: +SKIP
    Sum((Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*(Cs[i, j]*Cs[i, k]*dC_dTs[j, k] + Cs[i, j]*Cs[j, k]*dC_dTs[i, k] + Cs[i, k]*Cs[j, k]*dC_dTs[i, j])*KroneckerDelta(i, m)*zs[j]*zs[k]/(3*Cs[i, j]*Cs[i, k]*Cs[j, k]) + (Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*(Cs[i, j]*Cs[i, k]*dC_dTs[j, k] + Cs[i, j]*Cs[j, k]*dC_dTs[i, k] + Cs[i, k]*Cs[j, k]*dC_dTs[i, j])*KroneckerDelta(j, m)*zs[i]*zs[k]/(3*Cs[i, j]*Cs[i, k]*Cs[j, k]) + (Cs[i, j]*Cs[i, k]*Cs[j, k])**(1/3)*(Cs[i, j]*Cs[i, k]*dC_dTs[j, k] + Cs[i, j]*Cs[j, k]*dC_dTs[i, k] + Cs[i, k]*Cs[j, k]*dC_dTs[i, j])*KroneckerDelta(k, m)*zs[i]*zs[j]/(3*Cs[i, j]*Cs[i, k]*Cs[j, k]), (i, 0, nc), (j, 0, nc), (k, 0, nc))


    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    >>> zs = [.5, .3, .2]
    >>> d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs(zs, Cijs, dCij_dTs)
    [-1.5740994103e-11, -2.27267309501e-11, -3.56846953115e-11]
    '''
    N = len(zs)
    cC = [[0.0]*N for _ in range(N)] # numba: delete
#    cC = np.zeros((N, N)) # numba: uncomment

    for i in range(N):
        cC_row = cC[i]
        Ci_row = Cijs[i]
        for j in range(i):
            if Ci_row[j] > 0.0:
                cC_row[j] = cC[j][i] = Ci_row[j]**(1.0/3)/Ci_row[j]
            # else:
                # raise ValueError("Negative C")
        if Ci_row[i] > 0.0:
            cC_row[i] = Ci_row[i]**(1.0/3.0)/Ci_row[i]
        # else:
            # raise ValueError("Negative C")

    if d2C_dTdzs is None:
        d2C_dTdzs = [0.0]*N

    for i in range(N):
        i_sum = 0.0
        zi = zs[i]
        cCis = cC[i]
        for j in range(N):
            j_sum = 0.0
            zj = zs[j]
            zizj = zi*zj
            Cij = Cijs[i][j]
            cCij = cC[i][j]
            dCij = dCij_dTs[i][j]
            cCjs = cC[j]
            for k in range(N):
                t = cCij*cCis[k]*cCjs[k]
                c0 = (Cij*Cijs[i][k]*dCij_dTs[j][k]
                      + Cij*Cijs[j][k]*dCij_dTs[i][k]
                      + Cijs[i][k]*Cijs[j][k]*dCij)
                c1 = t*c0

                i_sum += c1*zj*zs[k]
                j_sum += c1*zi*zs[k]
                d2C_dTdzs[k] += c1*zizj

            d2C_dTdzs[j] += j_sum

        d2C_dTdzs[i] += i_sum

    for i in range(N):
        d2C_dTdzs[i] *= (1.0/3.0)
    return d2C_dTdzs

def CVirial_mixture_Orentlicher_Prausnitz(zs, Cijs):
    r'''Calculate the `C` third virial coefficient from a matrix of
    virial cross-coefficients. The diagonal is virial coefficients of the
    pure components.

    .. math::
        C = \sum_i \sum_j \sum_k y_i y_j y_k C_{ijk}(T)

    .. math::
        C_{ijk} = \left(C_{ij}C_{jk}C_{ik}\right)^{1/3}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]

    Returns
    -------
    C : float
        Third virial coefficient in density form [m^6/mol^2]

    Notes
    -----

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> zs = [.5, .3, .2]
    >>> CVirial_mixture_Orentlicher_Prausnitz(zs, Cijs)
    2.0790440095e-09

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    N = len(zs)
    Cij_cbrts = [[0.0]*N for _ in range(N)] # numba: delete
#    Cij_cbrts = np.zeros((N, N)) # numba: uncomment

    for i in range(N):
        Cij_cbrt_row = Cij_cbrts[i]
        Cij_row = Cijs[i]
        for j in range(i):
            if Cij_row[j] > 0.0:
                Cij_cbrt_row[j] = Cij_cbrts[j][i] = Cij_row[j]**(1.0/3)

        if Cij_row[i] > 0.0:
            Cij_cbrt_row[i] = Cij_row[i]**(1.0/3.0)

#     print(np.array(Cijs)**(1/3)/Cij_cbrts)
    C_diag, C_off = 0.0, 0.0
    # TODO: can we use symmetry to cut time down?
    for i in range(N):
        Cij_cbrts_i = Cij_cbrts[i]
        for j in range(N):
            x0 = zs[i]*zs[j]*Cij_cbrts_i[j]
            Cij_cbrts_j = Cij_cbrts[j]
            for k in range(j):
                if Cij_cbrts_i[k]*Cij_cbrts_j[k] > 0.0:
                    C_off += x0*zs[k]*Cij_cbrts_i[k]*Cij_cbrts_j[k]

            if Cij_cbrts_i[j]*Cij_cbrts_j[j] > 0.0:
                C_diag += x0*zs[j]*Cij_cbrts_i[j]*Cij_cbrts_j[j]
    return C_off*2.0 + C_diag

def dCVirial_mixture_dT_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs):
    r'''Calculate the first temperature derivative of the `C` third virial
    coefficient from matrices of
    virial cross-coefficients and their first temperature derivatives.

    .. math::
        \frac{\partial C}{\partial T} = \sum_i \sum_j \sum_k
        \frac{zi zj zk \sqrt[3]{\operatorname{Cij}{\left(T \right)} \operatorname{Cik}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)}} \left(\frac{\operatorname{Cij}{\left(T \right)}
        \operatorname{Cik}{\left(T \right)} \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}}{3}
        + \frac{\operatorname{Cij}{\left(T \right)} \operatorname{Cjk}{\left(T \right)} \frac{d}{d T}
        \operatorname{Cik}{\left(T \right)}}{3} + \frac{\operatorname{Cik}{\left(T \right)} \operatorname{Cjk}{\left(T \right)}
        \frac{d}{d T} \operatorname{Cij}{\left(T \right)}}{3}\right)}{\operatorname{Cij}{\left(T \right)}
        \operatorname{Cik}{\left(T \right)} \operatorname{Cjk}{\left(T \right)}}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    dCij_dTs : list[list[float]]
        First temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K]

    Returns
    -------
    dC_dT : float
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]

    Notes
    -----

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    >>> zs = [.5, .3, .2]
    >>> dCVirial_mixture_dT_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs)
    -7.2751517e-12

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    N = len(zs)
    Cij_pow_n23 = [[0.0]*N for _ in range(N)] # numba: delete
#    Cij_pow_n23 = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Cij_pow_n23_row = Cij_pow_n23[i]
        Cij_row = Cijs[i]
        for j in range(i):
            if Cij_row[j] > 0.0:
                Cij_pow_n23_row[j] = Cij_pow_n23[j][i] = Cij_row[j]**(1.0/3)/Cij_row[j]

        if Cij_row[i] > 0.0:
            Cij_pow_n23_row[i] = Cij_row[i]**(1.0/3.0)/Cij_row[i]

    # expect = np.array(Cijs)**(1/3)/Cijs
    # print(Cij_pow_n23/expect)

    dC = 0.0
    for i in range(N):
        for j in range(N):
            inner = 0.0
            x0 = zs[i]*zs[j]*Cij_pow_n23[i][j]
            if x0 < 0.0:
                continue
            Cij = Cijs[i][j]
            dCij = dCij_dTs[i][j]
            for k in range(N):
                t0 = Cij*Cijs[i][k]*dCij_dTs[j][k]
                t1 = Cij*Cijs[j][k]*dCij_dTs[i][k]
                t2 = Cijs[i][k]*Cijs[j][k]*dCij

                # Good
                # dC += zs[i]*zs[j]*zs[k]*(Cijs[i][j]*Cijs[i][k]*Cijs[j][k])**(1/3)*(t0 + t1 + t2)/(Cijs[i][j]*Cijs[i][k]*Cijs[j][k])
                # Good but suboptimal
                # dC += zs[i]*zs[j]*zs[k]*(Cijs[i][j]*Cijs[i][k]*Cijs[j][k])**(-2/3)*(t0 + t1 + t2)#/(Cijs[i][j]*Cijs[i][k]*Cijs[j][k])

                # Factor out the powers out
                term = Cij_pow_n23[i][k]*Cij_pow_n23[j][k]
                if term < 0.0:
                    continue
                inner += zs[k]*term*(t0 + t1 + t2)
            dC += inner*x0

    dC *= 1/3
    return dC

def d2CVirial_mixture_dT2_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s):
    r'''Calculate the second temperature derivative of the `C` third virial
    coefficient from matrices of
    virial cross-coefficients and their first and second temperature derivatives.

    .. math::
        \frac{\partial^2 C}{\partial T^2} = \sum_i \sum_j \sum_k z_i z_j z_k
        \frac{\sqrt[3]{\operatorname{Cij}{\left(T \right)} \operatorname{Cik}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)}} \left(\frac{\left(\operatorname{Cij}{\left(T \right)}
        \operatorname{Cik}{\left(T \right)} \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}
        + \operatorname{Cij}{\left(T \right)} \operatorname{Cjk}{\left(T \right)} \frac{d}{d T}
        \operatorname{Cik}{\left(T \right)} + \operatorname{Cik}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)} \frac{d}{d T} \operatorname{Cij}{\left(T \right)}\right)^{2}}
        {\operatorname{Cij}{\left(T \right)} \operatorname{Cik}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)}} - \frac{3 \left(\operatorname{Cij}{\left(T \right)}
        \operatorname{Cik}{\left(T \right)} \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}
        + \operatorname{Cij}{\left(T \right)} \operatorname{Cjk}{\left(T \right)} \frac{d}{d T}
        \operatorname{Cik}{\left(T \right)} + \operatorname{Cik}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)} \frac{d}{d T} \operatorname{Cij}{\left(T \right)}\right)
        \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}}{\operatorname{Cjk}{\left(T \right)}}
        - \frac{3 \left(\operatorname{Cij}{\left(T \right)} \operatorname{Cik}{\left(T \right)}
        \frac{d}{d T} \operatorname{Cjk}{\left(T \right)} + \operatorname{Cij}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)} \frac{d}{d T} \operatorname{Cik}{\left(T \right)}
        + \operatorname{Cik}{\left(T \right)} \operatorname{Cjk}{\left(T \right)} \frac{d}{d T}
        \operatorname{Cij}{\left(T \right)}\right) \frac{d}{d T} \operatorname{Cik}{\left(T \right)}}
        {\operatorname{Cik}{\left(T \right)}} - \frac{3 \left(\operatorname{Cij}{\left(T \right)}
        \operatorname{Cik}{\left(T \right)} \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}
        + \operatorname{Cij}{\left(T \right)} \operatorname{Cjk}{\left(T \right)} \frac{d}{d T}
        \operatorname{Cik}{\left(T \right)} + \operatorname{Cik}{\left(T \right)}
        \operatorname{Cjk}{\left(T \right)} \frac{d}{d T} \operatorname{Cij}{\left(T \right)}\right)
        \frac{d}{d T} \operatorname{Cij}{\left(T \right)}}{\operatorname{Cij}{\left(T \right)}}
        + 3 \operatorname{Cij}{\left(T \right)} \operatorname{Cik}{\left(T \right)} \frac{d^{2}}{d T^{2}}
        \operatorname{Cjk}{\left(T \right)} + 3 \operatorname{Cij}{\left(T \right)} \operatorname{Cjk}{\left(T \right)}
        \frac{d^{2}}{d T^{2}} \operatorname{Cik}{\left(T \right)} + 6 \operatorname{Cij}{\left(T \right)} \frac{d}{d T}
        \operatorname{Cik}{\left(T \right)} \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}
        + 3 \operatorname{Cik}{\left(T \right)} \operatorname{Cjk}{\left(T \right)}
        \frac{d^{2}}{d T^{2}} \operatorname{Cij}{\left(T \right)} + 6 \operatorname{Cik}{\left(T \right)}
        \frac{d}{d T} \operatorname{Cij}{\left(T \right)} \frac{d}{d T} \operatorname{Cjk}{\left(T \right)}
        + 6 \operatorname{Cjk}{\left(T \right)} \frac{d}{d T} \operatorname{Cij}{\left(T \right)}
        \frac{d}{d T} \operatorname{Cik}{\left(T \right)}\right)}{9 \operatorname{Cij}{\left(T \right)}
        \operatorname{Cik}{\left(T \right)} \operatorname{Cjk}{\left(T \right)}}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    dCij_dTs : list[list[float]]
        First temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K]
    d2Cij_dT2s : list[list[float]]
        Second temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K^2]

    Returns
    -------
    d2C_dT2 : float
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]

    Notes
    -----

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    >>> d2Cij_dT2s = [[ 2.6469e-14,  5.0512e-14,  1.1509e-13], [ 5.0512e-14,  9.3272e-14,  1.7836e-13], [ 1.1509e-13,  1.7836e-13, -1.4906e-13]]
    >>> zs = [.5, .3, .2]
    >>> d2CVirial_mixture_dT2_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s)
    6.7237107787e-14

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    N = len(zs)
    Cij_cbrts = [[0.0]*N for _ in range(N)] # numba: delete
#    Cij_cbrts = np.zeros((N, N)) # numba: uncomment
    # for i in range(N):
    #     for j in range(N):
    #         Cij_cbrts[i][j] = (Cijs[i][j])**(1/3)

    for i in range(N):
        Cij_cbrt_row = Cij_cbrts[i]
        Cij_row = Cijs[i]
        for j in range(i):
            if Cij_row[j] > 0.0:
                Cij_cbrt_row[j] = Cij_cbrts[j][i] = Cij_row[j]**(1.0/3)

        if Cij_row[i] > 0.0:
            Cij_cbrt_row[i] = Cij_row[i]**(1.0/3.0)

    d2C = 0.0
    for i in range(N):
        for j in range(N):
            t0 = zs[i]*zs[j]
            x8 = dCij_dTs[i][j]
            x0 = Cijs[i][j]
            x0_inv = 1.0/x0
            x16 = d2Cij_dT2s[i][j]
            cCij = Cij_cbrts[i][j]
            for k in range(N):
                x1 = Cijs[i][k]
                x2 = Cijs[j][k]
                x3 = x1*x2
                if x0*x3 < 0.0:
                    continue
                x4 = x0*x1
                x5 = x0*x2
                x6 = dCij_dTs[i][k]
                x7 = dCij_dTs[j][k]
                x9 = 6.0*x8
                x11 = x3*x8 + x4*x7 + x5*x6
                x12 = 3.0*x11
                x13 = 1.0/x1
                x14 = 1.0/x2
                x15 = x0_inv*x13*x14
                big = (x15*cCij*Cij_cbrts[i][k]*Cij_cbrts[j][k]*(6.0*x0*x6*x7 + x1*x7*x9 - x0_inv*x12*x8
                        + x11*x11*x15 - x12*x13*x6 - x12*x14*x7 + x2*x6*x9
                        + 3.0*x3*x16 + 3.0*x4*d2Cij_dT2s[j][k] + 3.0*x5*d2Cij_dT2s[i][k]))
                d2C += t0*zs[k]*big
    return d2C*(1/9.0)

def d3CVirial_mixture_dT3_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s,
                                                d3Cij_dT3s):
    r'''Calculate the third temperature derivative of the `C` third virial
    coefficient from matrices of
    virial cross-coefficients and their first, second, and third temperature
    derivatives.

    The expression is quite lengthy and not shown here [1]_.

    .. math::
        \frac{\partial^3 C}{\partial T^3}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of each species, [-]
    Cijs : list[list[float]]
        Third virial binary interaction coefficients in density form [m^6/mol^2]
    dCij_dTs : list[list[float]]
        First temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K]
    d2Cij_dT2s : list[list[float]]
        Second temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K^2]
    d3Cij_dT3s : list[list[float]]
        Third temperature derivative of third virial binary interaction
        coefficients in density form [m^6/mol^2/K^2^2]

    Returns
    -------
    d3C_dT3 : float
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Notes
    -----

    Examples
    --------
    >>> Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]
    >>> dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    >>> d2Cij_dT2s = [[ 2.6469e-14,  5.0512e-14,  1.1509e-13], [ 5.0512e-14,  9.3272e-14,  1.7836e-13], [ 1.1509e-13,  1.7836e-13, -1.4906e-13]]
    >>> d3Cij_dT3s = [[-4.2300e-16, -7.9727e-16, -1.6962e-15], [-7.9727e-16, -1.3826e-15, -1.4525e-15], [-1.6962e-15, -1.4525e-15,  1.9786e-14]]
    >>> zs = [.5, .3, .2]
    >>> d3CVirial_mixture_dT3_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s, d3Cij_dT3s)
    -3.7358368e-16

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    N = len(zs)
    d3C = 0.0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                x0 = Cijs[j][k]
                x1 = 1/x0
                x2 = Cijs[i][j]
                x3 = Cijs[i][k]
                x4 = x0*x3
                x2x4 = x2*x4

                # Create a discontinuity in this mixing rule if we were going to have complex components
                if x2x4 < 0.0:
                    continue
                x6 = d2Cij_dT2s[j][k]
                x7 = dCij_dTs[i][k]
                x9 = dCij_dTs[j][k]
                x10 = d2Cij_dT2s[i][k]
                x11 = dCij_dTs[i][j]
                x13 = d2Cij_dT2s[i][j]
                term = (x2x4)**(1/3)
                x5 = x2/3.0
                x8 = x2*x7
                x12 = x3*x6
                x14 = x3*x9
                x15 = x0*x10
                x16 = x0*x7
                x17 = 2.0*x11
                x18 = x7*x9
                x19 = 1.0/x2
                x20 = x11*x4 + x14*x2 + x16*x2
                x21 = x20/3.0
                x22 = x19*x21
                x23 = 1.0/x3
                x24 = x2**(-2.0)
                x25 = 2.0*x20/3.0
                x26 = x3**(-2.0)
                x27 = x0**(-2.0)
                x28 = x23*x7
                x29 = x11*x19
                x30 = x25*x29
                x31 = x1*x9
                x32 = x1*x23
                x33 = x20**2.0/3.0
                x34 = x1*x33
                x35 = x19*x23
                x36 = x12*x2 + x13*x4 + x14*x17 + x15*x2 + x16*x17 + 2.0*x8*x9
                x37 = 2.0*x36/3.0
                big = (x1*x35*term*(x0*x5*d3Cij_dT3s[i][k] - x1*x21*x6
                      + x10*x2*x9 - x10*x21*x23 + x11**2.0*x24*x25 + x11*x12 + x11*x15
                      - x11*x23*x24*x34 + x13*x14 + x13*x16 - x13*x22 + x17*x18 + x18*x25*x32
                      - x19*x26*x34*x7 + x20**3.0*x24*x26*x27/27.0 + x22*x32*x36 + x25*x26*x7**2.0
                      + x25*x27*x9**2.0 - x27*x33*x35*x9 + x28*x30 - x28*x37 - x29*x37
                      + x3*x5*d3Cij_dT3s[j][k] + x30*x31 - x31*x37 + x4*d3Cij_dT3s[i][j]/3.0 + x6*x8))
                d3C += zs[i]*zs[j]*zs[k]*big
    return d3C


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
        Third virial coefficient in density form [m^6/mol^2]
    dC_dT : float
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2 : float
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3 : float
        Third temperature derivative of third virial coefficient in density
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
    x1 = x0**(-21.0/2.0) # not worth optimizing, power tree is too large
    x2 = x0**(-14.0/5.0)
    x6 = -2000.0*x5
    x8 = 60.0*omega

    C = -x7*(2.0*omega*(114.0*x1 - 885.0*x2 + 150.0*x3 + x6 + 1338.0) + 313.0*x1 - 2432.0*x2 - 1407.0)*(1.0/100000.0)
    dC = x7*(32865.0*x1 - 68096.0*x2 + x8*(399.0*x1 - 826.0*x2 + 300.0*x3 + x6))*Tinv*(1.0/1000000.0)
    d2C = -x7*(3779475.0*x1 - 2587648.0*x2 + x8*(45885.0*x1 - 31388.0*x2 + 21000.0*x3 - 80000.0*x5))*Tinv2*(1.0/10000000.0)
    d3C = 3.0*x4*x7*(20.0*omega*(5735625.0*x1 - 1506624.0*x2 + 1680000.0*x3 - 4000000.0*x5) + 157478125.0*x1 - 41402368.0*x2)*(1.0/100000000.0)
    return C, dC, d2C, d3C

def CVirial_Orbey_Vera_vec(T, Tcs, Pcs, omegas, Cs=None, dC_dTs=None,
                         d2C_dT2s=None, d3C_dT3s=None):
    r'''Perform a vectorized calculation of the Orbey-Vera C virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Cs : list[float], optional
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[float], optional
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[float], optional
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[float], optional
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Returns
    -------
    Cs : list[float]
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[float]
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[float]
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[float]
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Cs is None:
        Cs = [0.0]*N
    if dC_dTs is None:
        dC_dTs = [0.0]*N
    if d2C_dT2s is None:
        d2C_dT2s = [0.0]*N
    if d3C_dT3s is None:
        d3C_dT3s = [0.0]*N
    for i in range(N):
        C, dC, d2C, d3C = CVirial_Orbey_Vera(T, Tcs[i], Pcs[i], omegas[i])
        Cs[i] = C
        dC_dTs[i] = dC
        d2C_dT2s[i] = d2C
        d3C_dT3s[i] = d3C
    return Cs, dC_dTs, d2C_dT2s, d3C_dT3s

def CVirial_Orbey_Vera_mat(T, Tcs, Pcs, omegas, Cs=None, dC_dTs=None,
                         d2C_dT2s=None, d3C_dT3s=None):
    r'''Perform a matrix calculation of the Orbey-Vera C virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Cs : list[list[float]], optional
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[list[float]], optional
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[list[float]], optional
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[list[float]], optional
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Returns
    -------
    Cs : list[list[float]]
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[list[float]]
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[list[float]]
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[list[float]]
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Cs is None:
        Cs = [[0.0]*N for _ in range(N)] # numba: delete
#        Cs = np.zeros((N, N)) # numba: uncomment
    if dC_dTs is None:
        dC_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dC_dTs = np.zeros((N, N)) # numba: uncomment
    if d2C_dT2s is None:
        d2C_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2C_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3C_dT3s is None:
        d3C_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3C_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        omega_row = omegas[i]

        C_row = Cs[i]
        dC_row = dC_dTs[i]
        d2C_row = d2C_dT2s[i]
        d3C_row = d3C_dT3s[i]

        for j in range(N):
            C, dC, d2C, d3C = CVirial_Orbey_Vera(T, Tc_row[j], Pc_row[j], omega_row[j])
            C_row[j] = C
            dC_row[j] = dC
            d2C_row[j] = d2C
            d3C_row[j] = d3C
    return Cs, dC_dTs, d2C_dT2s, d3C_dT3s

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
        Third virial coefficient in density form [m^6/mol^2]
    dC_dT : float
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2 : float
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3 : float
        Third temperature derivative of third virial coefficient in density
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


def CVirial_Liu_Xiang_vec(T, Tcs, Pcs, Vcs, omegas, Cs=None, dC_dTs=None,
                         d2C_dT2s=None, d3C_dT3s=None):
    r'''Perform a vectorized calculation of the Xiang C virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[float]
        Critical temperature of fluids [K]
    Pcs : list[float]
        Critical pressure of the fluids [Pa]
    Vcs : list[float]
        Critical volume of the fluids [m^3/mol]
    omegas : list[float]
        Acentric factor for fluids, [-]
    Cs : list[float], optional
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[float], optional
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[float], optional
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[float], optional
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Returns
    -------
    Cs : list[float]
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[float]
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[float]
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[float]
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Cs is None:
        Cs = [0.0]*N
    if dC_dTs is None:
        dC_dTs = [0.0]*N
    if d2C_dT2s is None:
        d2C_dT2s = [0.0]*N
    if d3C_dT3s is None:
        d3C_dT3s = [0.0]*N
    for i in range(N):
        C, dC, d2C, d3C = CVirial_Liu_Xiang(T, Tcs[i], Pcs[i], Vcs[i], omegas[i])
        Cs[i] = C
        dC_dTs[i] = dC
        d2C_dT2s[i] = d2C
        d3C_dT3s[i] = d3C
    return Cs, dC_dTs, d2C_dT2s, d3C_dT3s

def CVirial_Liu_Xiang_mat(T, Tcs, Pcs, Vcs, omegas, Cs=None, dC_dTs=None,
                         d2C_dT2s=None, d3C_dT3s=None):
    r'''Perform a matrix calculation of the Xiang C virial coefficient model
    and its first three temperature derivatives.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tcs : list[list[float]]
        Critical temperature of fluids [K]
    Pcs : list[list[float]]
        Critical pressure of the fluids [Pa]
    Vcs : list[list[float]]
        Critical volume of the fluids [m^3/mol]
    omegas : list[list[float]]
        Acentric factor for fluids, [-]
    Cs : list[list[float]], optional
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[list[float]], optional
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[list[float]], optional
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[list[float]], optional
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Returns
    -------
    Cs : list[list[float]]
        Third virial coefficient in density form [m^6/mol^2]
    dC_dTs : list[list[float]]
        First temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K]
    d2C_dT2s : list[list[float]]
        Second temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^2]
    d3C_dT3s : list[list[float]]
        Third temperature derivative of third virial coefficient in density
        form [m^6/mol^2/K^3]

    Notes
    -----
    '''
    N = len(Tcs)
    if Cs is None:
        Cs = [[0.0]*N for _ in range(N)] # numba: delete
#        Cs = np.zeros((N, N)) # numba: uncomment
    if dC_dTs is None:
        dC_dTs = [[0.0]*N for _ in range(N)] # numba: delete
#        dC_dTs = np.zeros((N, N)) # numba: uncomment
    if d2C_dT2s is None:
        d2C_dT2s = [[0.0]*N for _ in range(N)] # numba: delete
#        d2C_dT2s = np.zeros((N, N)) # numba: uncomment
    if d3C_dT3s is None:
        d3C_dT3s = [[0.0]*N for _ in range(N)] # numba: delete
#        d3C_dT3s = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        Tc_row = Tcs[i]
        Pc_row = Pcs[i]
        Vc_row = Vcs[i]
        omega_row = omegas[i]

        C_row = Cs[i]
        dC_row = dC_dTs[i]
        d2C_row = d2C_dT2s[i]
        d3C_row = d3C_dT3s[i]

        for j in range(N):
            C, dC, d2C, d3C = CVirial_Liu_Xiang(T, Tc_row[j], Pc_row[j], Vc_row[j], omega_row[j])
            C_row[j] = C
            dC_row[j] = dC
            d2C_row[j] = d2C
            d3C_row[j] = d3C
    return Cs, dC_dTs, d2C_dT2s, d3C_dT3s


### Mixing Rules

def Meng_Duan_2005_virial_CSP_kij_alkane(nci, ncj):
    # This whole thing can be cached as a lookup table
    # Also, not temperature dependent so no issues here
    if nci > ncj:
        nci, ncj = ncj, nci
    m = 0.00678/(1.0 + 0.336*nci)
    kij = m*(log(ncj - nci + 1.0)**(7.0/2.0))
    return kij

def Meng_Duan_2005_virial_CSP_kij_alkane_N2(nci):
    if nci < 1:
        return 0.0
    m = 0.04311
    return m*log(nci + 1.0)**1.5 # equation 15

def Meng_Duan_2005_virial_CSP_kij_alkane_CO2(nci):
    if nci < 1:
        return 0.0
    m = 0.07475
    return m*log(nci + 1.0)**1.5 # equation 15

CO2_CAS = '124-38-9'
N2_CAS = '7727-37-9'

def Meng_Duan_2005_virial_CSP_kijs(CASs, atomss):
    r'''Calculates a binary interaction parameter for the calculation of Bij
    binary virial coefficient as shown in [1]_. This implements a correlation
    of alkane-alkane, CO2-alkane, and N2-alkane.

    The equation this kij is used in is

    .. math::
        T_{cij} = \sqrt{T_{ci}T_{cj}}(1-k_{ij})

    Parameters
    ----------
    CASs : list[str]
        CAS registration numbers for each component, [-]
    atomss : list[dict]
        Breakdown of each component into its elements and their counts, as a
        dict, [-]

    Returns
    -------
    kijs : list[list[float]]
        Binary interaction parameters, [-]

    Notes
    -----

    Examples
    --------
    >>> CASs = ['74-82-8', '74-84-0', '124-38-9', '7727-37-9', '7439-89-6']
    >>> atomss = [{'C': 1, 'H': 4}, {'C': 2, 'H': 6}, {'C': 1, 'O': 2}, {'N': 2}, {'Fe': 1}]
    >>> kijs = Meng_Duan_2005_virial_CSP_kijs(CASs=CASs, atomss=atomss)

    References
    ----------
    .. [1] Meng, Long, and Yuan-Yuan Duan. "Prediction of the Second Cross
       Virial Coefficients of Nonpolar Binary Mixtures." Fluid Phase Equilibria
       238 (December 1, 2005): 229-38.
       https://doi.org/10.1016/j.fluid.2005.10.007.
    '''
    N = len(CASs)
    kijs = [[0.0]*N for _ in range(N)]
    for i in range(N):
        CAS1 = CASs[i]
        kij_row = kijs[i]
        C_i = atomss[i].get('C', 0)
        # symmetrical
        for j in range(i):
            CAS2 = CASs[j]

            if CAS1 == CO2_CAS:
                C = atomss[j].get('C', 0)
                kij_row[j] = kijs[j][i] = Meng_Duan_2005_virial_CSP_kij_alkane_CO2(C)
            elif CAS2 == CO2_CAS:
                kij_row[j] = kijs[j][i] = Meng_Duan_2005_virial_CSP_kij_alkane_CO2(C_i)

            elif CAS1 == N2_CAS:
                C = atomss[j].get('C', 0)
                kij_row[j] = kijs[j][i] = Meng_Duan_2005_virial_CSP_kij_alkane_N2(C)
            elif CAS2 == N2_CAS:
                kij_row[j] = kijs[j][i] = Meng_Duan_2005_virial_CSP_kij_alkane_N2(C_i)
            elif C_i and atomss[j].get('C', 0):
                kij_row[j] = kijs[j][i] = Meng_Duan_2005_virial_CSP_kij_alkane(C_i, atomss[j].get('C', 0))
            else:
                continue
    return kijs

def Tarakad_Danner_virial_CSP_kijs(Vcs):
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
    >>> Tarakad_Danner_virial_CSP_kijs(Vcs=[0.000168, 0.000316])
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
#    kijs = np.zeros((N, N)) # numba: uncomment
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
    >>> kijs = Tarakad_Danner_virial_CSP_kijs(Vcs=[0.000168, 0.000316])
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
#    Tcijs = np.zeros((N, N)) # numba: uncomment
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
    >>> kijs = Tarakad_Danner_virial_CSP_kijs(Vcs=[0.000168, 0.000316])
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
#    Pcijs = np.zeros((N, N)) # numba: uncomment
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
#    omegaijs = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        omegai = omegas[i]
        r = omegaijs[i]
        for j in range(N):
            r[j] = 0.5*(omegai + omegas[j])
    return omegaijs

def Lee_Kesler_virial_CSP_Vcijs(Vcs):
    r'''Calculates the corresponding states critical volumes for the
    calculation of Vcijs
    binary virial coefficient as shown in [1]_ and [2]_.

    .. math::
        V_{cij} = \frac{1}{8}\left(V_{c,i}^{1/3} + V_{c,j}^{1/3}
            \right)^3

    Parameters
    ----------
    Vcs : list[float]
        Critical volume of the fluids [m^3/mol]

    Returns
    -------
    Vcijs : list[list[float]]
        CSP critical volumes for each pair of species, [m^3/mol]

    Notes
    -----
    [1]_ cites this as Lee-Kesler rules.

    Examples
    --------
    >>> Lee_Kesler_virial_CSP_Vcijs(Vcs=[0.000168, 0.000316])
    [[0.000168, 0.00023426], [0.000234265, 0.000316]]

    References
    ----------
    .. [1] Estela-Uribe, J. F., and J. Jaramillo. "Generalised Virial Equation
       of State for Natural Gas Systems." Fluid Phase Equilibria 231, no. 1
       (April 1, 2005): 84-98. https://doi.org/10.1016/j.fluid.2005.01.005.
    .. [2] Lee, Byung Ik, and Michael G. Kesler. "A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States." AIChE
       Journal 21, no. 3 (1975): 510-27. https://doi.org/10.1002/aic.690210313.
    '''
    N = len(Vcs)
    Vcijs = [[0.0]*N for i in range(N)] # numba: delete
#    Vcijs = np.zeros((N, N)) # numba: uncomment
    Vc_cbrts = [0.0]*N
    for i in range(N):
        Vc_cbrts[i] = Vcs[i]**(1.0/3.0)
    for i in range(N):
        Vci_cbrt = Vc_cbrts[i]
        Vcij_row = Vcijs[i]

        for j in range(N):
            f = Vci_cbrt + Vc_cbrts[j]
            Vcij_row[j] = 0.125*f*f*f
    return Vcijs



def dV_dzs_virial(B, C, V, dB_dzs, dC_dzs, dV_dzs=None):
    r'''Calculates first mole fraction derivative of volume for the virial
    equation of state.

    .. math::
        \frac{\partial V}{\partial z_i} = \frac{V(V\frac{\partial B}{\partial z_i} + \frac{\partial C}{\partial z_i} )}{2BV + 3C + V^2}

    Parameters
    ----------
    B : float
        Second virial coefficient in density form [m^3/mol]
    C : float
        Third virial coefficient in density form [m^6/mol^2]
    V : float
        Molar volume from virial equation, [m^3/mol]
    dB_dzs : list[float]
        First mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]
    dC_dzs : list[float]
        First derivatives of C with respect to mole fraction, [m^6/mol^2]
    dV_dzs : list[float], optional
        Array for first derivatives of molar volume with respect to mole fraction, [m^3/mol]

    Returns
    -------
    dV_dzs : list[float]
        First derivatives of molar volume with respect to mole fraction, [m^3/mol]

    Notes
    -----
    This expression was derived with SymPy as follows:

    >>> from sympy import * # doctest: +SKIP
    >>> Z, R, T, P, z1 = symbols('Z, R, T, P, z1') # doctest: +SKIP
    >>> B, C, V = symbols('B, C, V', cls=Function) # doctest: +SKIP
    >>> base =Eq(P*V(z1)/(R*T), 1 + B(z1)/V(z1) + C(z1)/V(z1)**2) # doctest: +SKIP
    >>> P_sln = solve(base, P)[0] # doctest: +SKIP
    >>> solve(diff(P_sln, z1), Derivative(V(z1), z1)) # doctest: +SKIP
    [(V(z1)*Derivative(B(z1), z1) + Derivative(C(z1), z1))*V(z1)/(2*B(z1)*V(z1) + 3*C(z1) + V(z1)**2)]

    Examples
    --------
    >>> dV_dzs_virial(B=-5.130920247359858e-05, C=2.6627784284381213e-09, V=0.024892080086430797, dB_dzs=[-4.457911131778849e-05, -9.174964457681726e-05, -0.0001594258679841028], dC_dzs=[6.270599057032657e-09, 7.766612052069565e-09, 9.503031492910165e-09])
    [-4.4510120473455416e-05, -9.181495962913208e-05, -0.00015970040988493522]
    '''
    N = len(dB_dzs)
    if dV_dzs is None:
        dV_dzs = [0.0]*N
    for i in range(N):
        dV_dzs[i] = (V*dB_dzs[i] + dC_dzs[i])*V/(2.0*B*V + 3.0*C + V*V)
    return dV_dzs

def d2V_dzizjs_virial(B, C, V, dB_dzs, dC_dzs, dV_dzs, d2B_dzizjs, d2C_dzizjs,
                      d2V_dzizjs=None):
    r'''Calculates second mole fraction derivative of volume for the virial
    equation of state.

    .. math::
        \frac{\partial^2 V}{\partial z_i \partial z_j}

    Parameters
    ----------
    B : float
        Second virial coefficient in density form [m^3/mol]
    C : float
        Third virial coefficient in density form [m^6/mol^2]
    V : float
        Molar volume from virial equation, [m^3/mol]
    dB_dzs : list[float]
        First mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]
    dC_dzs : list[float]
        First derivatives of C with respect to mole fraction, [m^6/mol^2]
    dV_dzs : list[float]
        First derivatives of molar volume with respect to mole fraction, [m^3/mol]
    d2B_dzizjs : list[list[float]]
        Second mole fraction derivatives of second virial coefficient in
        density form [m^3/mol]
    d2C_dzizjs : list[list[float]]
        Second derivatives of C with respect to mole fraction, [m^6/mol^2]
    d2V_dzizjs : list[list[float]], optional
        Array for second derivatives of molar volume with respect to mole
        fraction, [m^3/mol]

    Returns
    -------
    d2V_dzizjs : list[list[float]]
        Second derivatives of molar volume with respect to mole
        fraction, [m^3/mol]

    Notes
    -----
    This expression was derived with SymPy as follows:

    >>> from sympy import * # doctest: +SKIP
    >>> Z, R, T, P, z1 = symbols('Z, R, T, P, z1') # doctest: +SKIP
    >>> B, C, V = symbols('B, C, V', cls=Function) # doctest: +SKIP
    >>> base =Eq(P*V(z1)/(R*T), 1 + B(z1)/V(z1) + C(z1)/V(z1)**2) # doctest: +SKIP
    >>> P_sln = solve(base, P)[0] # doctest: +SKIP
    >>> solve(diff(P_sln, z1), Derivative(V(z1), z1)) # doctest: +SKIP
    [(V(z1)*Derivative(B(z1), z1) + Derivative(C(z1), z1))*V(z1)/(2*B(z1)*V(z1) + 3*C(z1) + V(z1)**2)]

    Examples
    --------
    >>> d2C_dzizjs = [[1.0287075724127612e-08, 1.2388277824773021e-08, 1.4298813522844275e-08], [1.2388277824773021e-08, 1.514162073913238e-08, 1.8282527232061114e-08], [1.4298813522844275e-08, 1.8282527232061114e-08, 2.3350122217403063e-08]]
    >>> d2B_dzizjs = [[-1.0639357784985337e-05, -3.966321845899801e-05, -7.53987684376414e-05], [-3.966321845899801e-05, -8.286257232134107e-05, -0.00014128571574782375], [-7.53987684376414e-05, -0.00014128571574782375, -0.00024567752140887547]]
    >>> dB_dzs = [-4.457911131778849e-05, -9.174964457681726e-05, -0.0001594258679841028]
    >>> dC_dzs = [6.270599057032657e-09, 7.766612052069565e-09, 9.503031492910165e-09]
    >>> dV_dzs = [-4.4510120473455416e-05, -9.181495962913208e-05, -0.00015970040988493522]
    >>> d2V_dzizjs_virial(B=-5.130920247359858e-05, C=2.6627784284381213e-09, V=0.024892080086430797, dB_dzs=dB_dzs, dC_dzs=dC_dzs, dV_dzs=dV_dzs, d2B_dzizjs=d2B_dzizjs, d2C_dzizjs=d2C_dzizjs)
    [[-1.04268917389e-05, -3.9654694588e-05, -7.570310078e-05], [-3.9654694588e-05, -8.3270116767e-05, -0.0001423083584], [-7.5703100789e-05, -0.000142308358, -0.00024779788]]
    '''
    N = len(dB_dzs)
    if d2V_dzizjs is None:
        d2V_dzizjs = [[0.0]*N for _ in range(N)] # numba: delete
        # d2V_dzizjs = np.zeros((N, N)) # numba: uncomment
    for i in range(N):
        for j in range(N):
            x0 = V
            x1 = C
            x2 = x0*x0
            x3 = B
            x4 = x0*x3
            x5 = dV_dzs[j]
            x6 = dV_dzs[i]
            x7 = x5*x6
            x8 = 3*x0
            x9 = 2*x2
            x10 = x5*x9
            d2V_dzizjs[i][j] = (x0*x0*x0*d2B_dzizjs[i][j]+ 12.0*x1*x7 + x10*x6 - x10*dB_dzs[i]
                                + x2*d2C_dzizjs[i][j] + 6.0*x4*x7 - x5*x8*dC_dzs[i] - x6*x8*dC_dzs[j]
                                - x6*x9*dB_dzs[j])/(x0*(3.0*x1 + x2 + 2.0*x4))
    return d2V_dzizjs
