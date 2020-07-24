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
from math import exp, log

from fluids.numerics import newton, brenth, oscillation_checker, secant, NotBoundedError

from chemicals.rachford_rice import flash_inner_loop

__all__ = ['K_value','Wilson_K_value', 'PR_water_K_value', 'flash_wilson', 
           'flash_Tb_Tc_Pc', 'flash_ideal']

__numba_additional_funcs__ = ['err_Wilson_TVF', 'err_Wilson_PVF']


def K_value(P=None, Psat=None, phi_l=None, phi_g=None, gamma=None, Poynting=1.0):
    r'''Calculates the equilibrium K-value assuming Raoult's law,
    or an equation of state model, or an activity coefficient model,
    or a combined equation of state-activity model.

    The calculation procedure will use the most advanced approach with the
    provided inputs:

        * If `P`, `Psat`, `phi_l`, `phi_g`, and `gamma` are provided, use the
          combined approach.
        * If `P`, `Psat`, and `gamma` are provided, use the modified Raoult's
          law.
        * If `phi_l` and `phi_g` are provided, use the EOS only method.
        * If `P` and `Psat` are provided, use Raoult's law.

    Definitions:

    .. math::
        K_i=\frac{y_i}{x_i}

    Raoult's law:

    .. math::
        K_i = \frac{P_{i}^{sat}}{P}

    Activity coefficient, no EOS (modified Raoult's law):

    .. math::
        K_i = \frac{\gamma_i P_{i}^{sat}}{P}

    Equation of state only:

    .. math::
        K_i = \frac{\phi_i^l}{\phi_i^v} = \frac{f_i^l y_i}{f_i^v x_i}

    Combined approach (liquid reference fugacity coefficient is normally
    calculated the saturation pressure for it as a pure species; vapor fugacity
    coefficient calculated normally):

    .. math::
        K_i = \frac{\gamma_i P_i^{sat} \phi_i^{l,ref}}{\phi_i^v P}

    Combined approach, with Poynting Correction Factor (liquid molar volume in
    the integral is for i as a pure species only):

    .. math::
        K_i = \frac{\gamma_i P_i^{sat} \phi_i^{l, ref} \exp\left[\frac{
        \int_{P_i^{sat}}^P V_i^l dP}{RT}\right]}{\phi_i^v P}

    Parameters
    ----------
    P : float
        System pressure, optional
    Psat : float
        Vapor pressure of species i, [Pa]
    phi_l : float
        Fugacity coefficient of species i in the liquid phase, either
        at the system conditions (EOS-only case) or at the saturation pressure
        of species i as a pure species (reference condition for the combined
        approach), optional [-]
    phi_g : float
        Fugacity coefficient of species i in the vapor phase at the system
        conditions, optional [-]
    gamma : float
        Activity coefficient of species i in the liquid phase, optional [-]
    Poynting : float
        Poynting correction factor, optional [-]

    Returns
    -------
    K : float
        Equilibrium K value of component i, calculated with an approach
        depending on the provided inputs [-]

    Notes
    -----
    The Poynting correction factor is normally simplified as follows, due to
    a liquid's low pressure dependency:

    .. math::
        K_i = \frac{\gamma_i P_i^{sat} \phi_i^{l, ref} \exp\left[\frac{V_l
        (P-P_i^{sat})}{RT}\right]}{\phi_i^v P}

    Examples
    --------
    Raoult's law:

    >>> K_value(101325, 3000.)
    0.029607698001480384

    Modified Raoult's law:

    >>> K_value(P=101325, Psat=3000, gamma=0.9)
    0.026646928201332347

    EOS-only approach:

    >>> K_value(phi_l=1.6356, phi_g=0.88427)
    1.8496613025433408

    Gamma-phi combined approach:

    >>> K_value(P=1E6, Psat=1938800, phi_l=1.4356, phi_g=0.88427, gamma=0.92)
    2.8958055544121137

    Gamma-phi combined approach with a Poynting factor:

    >>> K_value(P=1E6, Psat=1938800, phi_l=1.4356, phi_g=0.88427, gamma=0.92,
    ... Poynting=0.999)
    2.8929097488577016

    References
    ----------
    .. [1] Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey.
       Chemical Thermodynamics for Process Simulation. 1st edition. Weinheim:
       Wiley-VCH, 2012.
    .. [2] Skogestad, Sigurd. Chemical and Energy Process Engineering. 1st
       edition. Boca Raton, FL: CRC Press, 2008.
    '''
    try:
        if gamma is not None:
            if phi_l is not None:
                return gamma*Psat*phi_l*Poynting/(phi_g*P)
            return gamma*Psat*Poynting/P
        elif phi_l is not None:
            return phi_l/phi_g
        return Psat/P
    except:
        raise TypeError('Input must consist of one set from (P, Psat, phi_l, '
                        'phi_g, gamma), (P, Psat, gamma), (phi_l, phi_g), (P, Psat)')


def Wilson_K_value(T, P, Tc, Pc, omega):
    r'''Calculates the equilibrium K-value for a component using Wilson's
    heuristic mode. This is very useful for initialization of stability tests
    and flashes.

    .. math::
        K_i = \frac{P_c}{P} \exp\left(5.37(1+\omega)\left[1 - \frac{T_c}{T}
        \right]\right)

    Parameters
    ----------
    T : float
        System temperature, [K]
    P : float
        System pressure, [Pa]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    K : float
        Equilibrium K value of component, calculated via the Wilson heuristic
        [-]

    Notes
    -----
    There has been little literature exploration of other formlulas for the
    same purpose. This model may be useful even for activity coefficient
    models.

    Note the K-values are independent of composition; the correlation is
    applicable up to 3.5 MPa.

    A description for how this function was generated can be found in [2]_.

    Examples
    --------
    Ethane at 270 K and 76 bar:

    >>> Wilson_K_value(270.0, 7600000.0, 305.4, 4880000.0, 0.098)
    0.2963932297479371

    The "vapor pressure" predicted by this equation can be calculated by
    multiplying by pressure:

    >>> Wilson_K_value(270.0, 7600000.0, 305.4, 4880000.0, 0.098)*7600000.0
    2252588.546084322

    References
    ----------
    .. [1] Wilson, Grant M. "A Modified Redlich-Kwong Equation of State,
       Application to General Physical Data Calculations." In 65th National
       AIChE Meeting, Cleveland, OH, 1969.
    .. [2] Peng, Ding-Yu, and Donald B. Robinson. "Two and Three Phase
       Equilibrium Calculations for Systems Containing Water." The Canadian
       Journal of Chemical Engineering, December 1, 1976.
       https://doi.org/10.1002/cjce.5450540620.
    '''
    return Pc/P*exp((5.37*(1.0 + omega)*(1.0 - Tc/T)))


def PR_water_K_value(T, P, Tc, Pc):
    r'''Calculates the equilibrium K-value for a component against water
    according to the Peng and Robinson (1976) heuristic.

    .. math::
        K_i = 10^6 \frac{P_{ri}}{T_{ri}}

    Parameters
    ----------
    T : float
        System temperature, [K]
    P : float
        System pressure, [Pa]
    Tc : float
        Critical temperature of chemical [K]
    Pc : float
        Critical pressure of chemical [Pa]

    Returns
    -------
    K : float
        Equilibrium K value of component with water as the other phase (
        not as the reference), calculated via this heuristic [-]

    Notes
    -----
    Note the K-values are independent of composition.

    Examples
    --------
    Octane at 300 K and 1 bar:

    >>> PR_water_K_value(300, 1e5, 568.7, 2490000.0)
    76131.19143239626

    References
    ----------
    .. [1] Peng, Ding-Yu, and Donald B. Robinson. "Two and Three Phase
       Equilibrium Calculations for Systems Containing Water." The Canadian
       Journal of Chemical Engineering, December 1, 1976.
       https://doi.org/10.1002/cjce.5450540620.
    '''
    Tr = T/Tc
    Pr = P/Pc
    return 1e6*Pr/Tr

def err_Wilson_TVF(P, N, VF, zs, K_Ps):
    P_inv = 1.0/P
    err, derr = 0.0, 0.0
    for i in range(N):
        x50 = K_Ps[i]*P_inv
        x0 = x50 - 1.0
        x1 = VF*x0
        x2 = 1.0/(x1 + 1.0)
        x3 = x2*zs[i]
        err += x0*x3
        derr += x50*P_inv*x3*(x1*x2 - 1.0)
    return err, derr

def err_Wilson_PVF(T_guess, N, P_inv, VF, Tcs, Pcs, Ks, zs, xs, x50s):
    err, derr = 0.0, 0.0
    T_inv = 1.0/T_guess
    T_inv2 = T_inv*T_inv
    for i in range(N):
        Ks[i] = Pcs[i]*exp(x50s[i]*(1.0 - Tcs[i]*T_inv))*P_inv
        dKi_dT = Ks[i]*x50s[i]*T_inv2*Tcs[i]
        x1 = Ks[i] - 1.0
        x2 = VF*x1
        x3 = 1.0/(x2 + 1.0)
        xs[i] = x3*zs[i]
        err += x1*xs[i]
        derr += xs[i]*(1.0 - x2*x3)*dKi_dT
    return err, derr

def flash_wilson(zs, Tcs, Pcs, omegas, T=None, P=None, VF=None):
    r'''PVT flash model using Wilson's equation - useful for obtaining initial
    guesses for more rigorous models, or it can be used as its own model.
    Capable of solving with two of `T`, `P`, and `VF` for the other one;
    that results in three solve modes, but for `VF=1` and `VF=0`, there are
    additional solvers; for a total of seven solvers implemented.

    This model uses `flash_inner_loop` to solve the Rachford-Rice problem.

    .. math::
        K_i = \frac{P_c}{P} \exp\left(5.37(1+\omega)\left[1 - \frac{T_c}{T}
        \right]\right)

    Parameters
    ----------
    zs : list[float]
        Mole fractions of the phase being flashed, [-]
    Tcs : list[float]
        Critical temperatures of all species, [K]
    Pcs : list[float]
        Critical pressures of all species, [Pa]
    omegas : list[float]
        Acentric factors of all species, [-]
    T : float, optional
        Temperature, [K]
    P : float, optional
        Pressure, [Pa]
    VF : float, optional
        Molar vapor fraction, [-]

    Returns
    -------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    VF : float
        Molar vapor fraction, [-]
    xs : list[float]
        Mole fractions of liquid phase, [-]
    ys : list[float]
        Mole fractions of vapor phase, [-]

    Notes
    -----
    For the cases where `VF` is 1 or 0 and T is known, an explicit solution is
    used. For the same cases where `P` and `VF` are known, there is no explicit
    solution available.

    There is an internal `Tmax` parameter, set to 50000 K; which, in the event
    of convergence of the Secant method, is used as a bounded for a bounded
    solver. It is used in the PVF solvers. This typically allows pressures
    up to 2 GPa to be converged to. However, for narrow-boiling mixtures, the
    PVF failure may occur at much lower pressures.

    Examples
    --------
    >>> Tcs = [305.322, 540.13]
    >>> Pcs = [4872200.0, 2736000.0]
    >>> omegas = [0.099, 0.349]
    >>> zs = [0.4, 0.6]
    >>> flash_wilson(zs=zs, Tcs=Tcs, Pcs=Pcs, omegas=omegas, T=300, P=1e5)
    (300, 100000.0, 0.42219453293637355, [0.020938815080034565, 0.9790611849199654], [0.9187741856225791, 0.08122581437742094])
    '''
    T_MAX = 50000.0
    N = len(zs)
    # Assume T and P to begin with
    if T is not None and P is not None:
        # numba is failing its type inferences
        P_inv = 1.0/P
        T_inv = 1.0/T
        Ks = [0.0]*N
        for i in range(N):
            Ks[i] = P_inv*Pcs[i]*exp((5.37*(1.0 + omegas[i])*(1.0 - Tcs[i]*T_inv)))
#        all_under_1, all_over_1 = True, True
#        for K in Ks:
#            if K < 1.0:
#                all_over_1 = False
#            else:
#                all_under_1 = False
#            if all_over_1:
#                raise ValueError("Fail")
#            elif all_under_1:
#                raise ValueError("Fail")
        ans = (T, P) + flash_inner_loop(zs=zs, Ks=Ks)
        return ans
    elif T is not None and VF is not None and VF == 0.0:
        ys = [0.0]*N
        P_bubble = 0.0
        T_inv = 1.0/T
        for i in range(N):
            v = zs[i]*Pcs[i]*exp((5.37*(1.0 + omegas[i])*(1.0 - Tcs[i]*T_inv)))
            P_bubble += v
            ys[i] = v
        P_inv = 1.0/P_bubble
        for i in range(N):
            ys[i] *= P_inv
        return (T, P_bubble, 0.0, zs, ys)
    elif T is not None and VF is not None and VF == 1.0:
        xs = [0.0]*N
        P_dew = 0.
        T_inv = 1.0/T
        for i in range(N):
            v = zs[i]/(Pcs[i]*exp((5.37*(1.0 + omegas[i])*(1.0 - Tcs[i]*T_inv))))
            P_dew += v
            xs[i] = v
        P_dew = 1./P_dew
        for i in range(N):
            xs[i] *= P_dew
        return (T, P_dew, 1.0, xs, zs)
    elif T is not None and VF is not None:
        # Solve for the pressure to create the desired vapor fraction
        P_bubble = 0.0
        P_dew = 0.
        T_inv = 1.0/T
        K_Ps = [0.0]*N
        for i in range(N):
            K_P = Pcs[i]*exp((5.37*(1.0 + omegas[i])*(1.0 - Tcs[i]*T_inv)))
            P_bubble += zs[i]*K_P
            P_dew += zs[i]/K_P
            K_Ps[i] = K_P
        P_dew = 1./P_dew
        '''Rachford-Rice esque solution in terms of pressure.
        from sympy import *
        N = 1
        cmps = range(N)
        zs = z0, z1, z2, z3 = symbols('z0, z1, z2, z3')
        Ks_P = K0_P, K1_P, K2_P, K3_P = symbols('K0_P, K1_P, K2_P, K3_P')
        VF, P = symbols('VF, P')
        tot = 0
        for i in cmps:
            tot += zs[i]*(Ks_P[i]/P - 1)/(1 + VF*(Ks_P[i]/P - 1))
        cse([tot, diff(tot, P)], optimizations='basic')
        '''
        P_guess = P_bubble + VF*(P_dew - P_bubble) # Linear interpolation
        P_calc = newton(err_Wilson_TVF, P_guess, fprime=True, bisection=True,
                   low=P_dew, high=P_bubble, args=(N, VF, zs, K_Ps))
        P_inv = 1.0/P_calc
        xs = K_Ps
        ys = [0.0]*N
        for i in range(N):
            Ki = K_Ps[i]*P_inv
            xi = zs[i]/(1.0 + VF*(Ki - 1.0))
            ys[i] = Ki*xi
            xs[i] = xi
        return (T, P_calc, VF, xs, ys)
    elif P is not None and VF is not None:
        P_inv = 1.0/P
        Ks = [0.0]*N
        xs = [0.0]*N
        x50s = [5.37]*N
        for i in range(N):
            x50s[i] *= omegas[i] + 1.0
        T_low, T_high = 1e100, 0.0
        logP = log(P)
        for i in range(N):
            T_K_1 = Tcs[i]*x50s[i]/(x50s[i] - logP + log(Pcs[i]))
            if T_K_1 < T_low:
                T_low = T_K_1
            if T_K_1 > T_high:
                T_high = T_K_1
        if T_low < 0.0:
            T_low = 1e-12
        if T_high <= 0.0:
            raise ValueError("No temperature exists which makes Wilson K factor above 1 - decrease pressure")
        if T_high < 0.1*T_MAX:
            T_guess = 0.5*(T_low + T_high)
        else:
            T_guess = 0.0
            for i in range(N):
                T_guess += zs[i]*Tcs[i]
            T_guess *= 0.666666
            if T_guess < T_low:
                T_guess = T_low + 1.0 # Take a nominal step
        T_calc = newton(err_Wilson_PVF, T_guess, fprime=True, low=T_low, xtol=1e-13, bisection=True, args=(N, P_inv, VF, Tcs, Pcs, Ks, zs, xs, x50s)) # High bound not actually a bound, only low bound
        if 1e-10 < T_calc < T_MAX:
            ys = x50s
            for i in range(N):
                ys[i] = xs[i]*Ks[i]
            return (T_calc, P, VF, xs, ys)

#    # Old code - may converge where the other will not
#    if P is not None and VF == 1.0:
#        def to_solve(T_guess):
#            # Avoid some nasty unpleasantness in newton
#            T_guess = abs(T_guess)
#            P_dew = 0.
#            for i in range(len(zs)):
#                P_dew += zs[i]/(Pcs[i]*exp((5.37*(1.0 + omegas[i])*(1.0 - Tcs[i]/T_guess))))
#            P_dew = 1./P_dew
##            print(P_dew - P, T_guess)
#            return P_dew - P
#        # 2/3 average critical point
#        T_guess = sum([.666*Tcs[i]*zs[i] for i in cmps])
#        try:
#            T_dew = abs(newton(to_solve, T_guess, maxiter=50, ytol=1e-2))
#        except Exception as e:
##            print(e)
#            T_dew = None
#        if T_dew is None or T_dew > T_MAX*5.0:
#            # Went insanely high T, bound it with brenth
#            T_low_guess = sum([.1*Tcs[i]*zs[i] for i in cmps])
#            try:
#                T_dew = brenth(to_solve, T_MAX, T_low_guess)
#            except NotBoundedError:
#                raise Exception("Bisecting solver could not find a solution between %g K and %g K" %(T_MAX, T_low_guess))
#        return flash_wilson(zs, Tcs, Pcs, omegas, T=T_dew, P=P)
#    elif P is not None and VF == 0.0:
#        def to_solve(T_guess):
#            T_guess = abs(T_guess)
#            P_bubble = 0.0
#            for i in cmps:
#                P_bubble += zs[i]*Pcs[i]*exp((5.37*(1.0 + omegas[i])*(1.0 - Tcs[i]/T_guess)))
#            return P_bubble - P
#        # 2/3 average critical point
#        T_guess = sum([.55*Tcs[i]*zs[i] for i in cmps])
#        try:
#            T_bubble = abs(newton(to_solve, T_guess, maxiter=50, ytol=1e-2))
#        except Exception as e:
#            T_bubble = None
#        if T_bubble is None or T_bubble > T_MAX*5.0:
#            # Went insanely high T, bound it with brenth
#            T_low_guess = sum([.1*Tcs[i]*zs[i] for i in cmps])
#            try:
#                T_bubble = brenth(to_solve, T_MAX, T_low_guess)
#            except NotBoundedError:
#                raise Exception("Bisecting solver could not find a solution between %g K and %g K" %(T_MAX, T_low_guess))
#
#        return flash_wilson(zs, Tcs, Pcs, omegas, T=T_bubble, P=P)
#    elif P is not None and VF is not None:
#        # Solve for in the middle of Pdew
#        T_low = flash_wilson(zs, Tcs, Pcs, omegas, P=P, VF=1)[0]
#        T_high = flash_wilson(zs, Tcs, Pcs, omegas, P=P, VF=0)[0]
##        print(T_low, T_high)
#        info = []
#        def err(T):
#            T_calc, P_calc, VF_calc, xs, ys = flash_wilson(zs, Tcs, Pcs, omegas, T=T, P=P)
##            if abs(VF_calc) > 100: # Did not work at all
##                VF_calc = abs(VF_calc)
#            info[:] = T_calc, P_calc, VF_calc, xs, ys
##            print(T, VF_calc - VF)
#            return VF_calc - VF
#        # Nasty function for tolerance; the default works and is good enough, could remove some
#        # iterations in the fuure
#        P = brenth(err, T_low, T_high, xtol=1e-14)
#        return tuple(info)
    else:
        raise ValueError("Provide two of P, T, and VF")


def flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=None, P=None, VF=None):
    r'''PVT flash model using a model published in [1]_, which provides a PT
    surface using only each compound's boiling temperature and critical
    temperature and pressure. This is useful for obtaining initial
    guesses for more rigorous models, or it can be used as its own model.
    Capable of solving with two of `T`, `P`, and `VF` for the other one;
    that results in three solve modes, but for `VF=1` and `VF=0`, there are
    additional solvers; for a total of seven solvers implemented.

    This model uses `flash_inner_loop` to solve the Rachford-Rice problem.

    .. math::
        K_i = \frac{P_{c,i}^{\left(\frac{1}{T} - \frac{1}{T_{b,i}} \right) /
        \left(\frac{1}{T_{c,i}} - \frac{1}{T_{b,i}} \right)}}{P}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of the phase being flashed, [-]
    Tbs : list[float]
        Boiling temperatures of all species, [K]
    Tcs : list[float]
        Critical temperatures of all species, [K]
    Pcs : list[float]
        Critical pressures of all species, [Pa]
    T : float, optional
        Temperature, [K]
    P : float, optional
        Pressure, [Pa]
    VF : float, optional
        Molar vapor fraction, [-]

    Returns
    -------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    VF : float
        Molar vapor fraction, [-]
    xs : list[float]
        Mole fractions of liquid phase, [-]
    ys : list[float]
        Mole fractions of vapor phase, [-]

    Notes
    -----
    For the cases where `VF` is 1 or 0 and T is known, an explicit solution is
    used. For the same cases where `P` and `VF` are known, there is no explicit
    solution available.

    There is an internal `Tmax` parameter, set to 50000 K; which, in the event
    of convergence of the Secant method, is used as a bounded for a bounded
    solver. It is used in the PVF solvers. This typically allows pressures
    up to 2 MPa to be converged to. Failures may still occur for other
    conditions.
    
    This model is based on [1]_, which aims to estimate dew and bubble points
    using the same K value formulation as used here. While this implementation
    uses a numerical solver to provide an exact bubble/dew point estimate,
    [1]_ suggests a sequential substitution and flowchart based solver with 
    loose tolerances. That model was also implemented, but found to be slower
    and less reliable than this implementation.
    

    Examples
    --------
    >>> Tcs = [305.322, 540.13]
    >>> Pcs = [4872200.0, 2736000.0]
    >>> Tbs = [184.55, 371.53]
    >>> zs = [0.4, 0.6]
    >>> flash_Tb_Tc_Pc(zs=zs, Tcs=Tcs, Pcs=Pcs, Tbs=Tbs, T=300, P=1e5)
    (300, 100000.0, 0.3807040748145384, [0.031157843036568357, 0.9688421569634317], [0.9999999998827085, 1.1729141887515062e-10])

    References
    ----------
    .. [1] Kandula, Vamshi Krishna, John C. Telotte, and F. Carl Knopf. "It’s 
       Not as Easy as It Looks: Revisiting Peng—Robinson Equation of State 
       Convergence Issues for Dew Point, Bubble Point and Flash Calculations."
       International Journal of Mechanical Engineering Education 41, no. 3 
       (July 1, 2013): 188-202. https://doi.org/10.7227/IJMEE.41.3.2.
    '''
    T_MAX = 50000
    N = len(zs)
    cmps = range(N)
    # Assume T and P to begin with
    if T is not None and P is not None:
        Ks = [Pcs[i]**((1.0/T - 1.0/Tbs[i])/(1.0/Tcs[i] - 1.0/Tbs[i]))/P for i in cmps]
        return (T, P) + flash_inner_loop(zs=zs, Ks=Ks, check=True)


    if T is not None and VF == 0:
        P_bubble = 0.0
        for i in cmps:
            P_bubble += zs[i]*Pcs[i]**((1.0/T - 1.0/Tbs[i])/(1.0/Tcs[i] - 1.0/Tbs[i]))
        return flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T, P=P_bubble)
    if T is not None and VF == 1:
        # Checked to be working vs. PT implementation.
        P_dew = 0.
        for i in cmps:
            P_dew += zs[i]/( Pcs[i]**((1.0/T - 1.0/Tbs[i])/(1.0/Tcs[i] - 1.0/Tbs[i])) )
        P_dew = 1./P_dew
        return flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T, P=P_dew)
    elif T is not None and VF is not None:
        # Solve for in the middle of Pdew
        P_low = flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T, VF=1)[1]
        P_high = flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T, VF=0)[1]
        info = []
        def err(P):
            T_calc, P_calc, VF_calc, xs, ys = flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T, P=P)
            info[:] = T_calc, P_calc, VF_calc, xs, ys
            return VF_calc - VF
        P = brenth(err, P_low, P_high)
        return tuple(info)

    elif P is not None and VF == 1:
        checker = oscillation_checker()
        def to_solve(T_guess):
            T_guess = abs(T_guess)
            P_dew = 0.
            for i in range(len(zs)):
                P_dew += zs[i]/( Pcs[i]**((1.0/T_guess - 1.0/Tbs[i])/(1.0/Tcs[i] - 1.0/Tbs[i])) )
            P_dew = 1./P_dew
            err = P_dew - P
            if checker(T_guess, err):
                raise ValueError("Oscillation")
#            print(T_guess, err)
            return err

        Tc_pseudo = sum([Tcs[i]*zs[i] for i in cmps])
        T_guess = 0.666*Tc_pseudo
        try:
            T_dew = abs(secant(to_solve, T_guess, maxiter=50, ytol=1e-2)) # , high=Tc_pseudo*3
        except:
            T_dew = None
        if T_dew is None or T_dew > T_MAX*5.0:
            # Went insanely high T, bound it with brenth
            T_low_guess = sum([.1*Tcs[i]*zs[i] for i in cmps])
            checker = oscillation_checker(both_sides=True, minimum_progress=.05)
            try:
                T_dew = brenth(to_solve, T_MAX, T_low_guess)
            except NotBoundedError:
                raise Exception("Bisecting solver could not find a solution between %g K and %g K" %(T_MAX, T_low_guess))
        return flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T_dew, P=P)

    elif P is not None and VF == 0:
        checker = oscillation_checker()
        def to_solve(T_guess):
            T_guess = abs(T_guess)
            P_bubble = 0.0
            for i in cmps:
                P_bubble += zs[i]*Pcs[i]**((1.0/T_guess - 1.0/Tbs[i])/(1.0/Tcs[i] - 1.0/Tbs[i]))

            err = P_bubble - P
            if checker(T_guess, err):
                raise ValueError("Oscillation")

#            print(T_guess, err)
            return err
        # 2/3 average critical point
        Tc_pseudo = sum([Tcs[i]*zs[i] for i in cmps])
        T_guess = 0.55*Tc_pseudo
        try:
            T_bubble = abs(secant(to_solve, T_guess, maxiter=50, ytol=1e-2)) # , high=Tc_pseudo*4
        except Exception as e:
#            print(e)
            checker = oscillation_checker(both_sides=True, minimum_progress=.05)
            T_bubble = None
        if T_bubble is None or T_bubble > T_MAX*5.0:
            # Went insanely high T (or could not converge because went too high), bound it with brenth
            T_low_guess = 0.1*Tc_pseudo
            try:
                T_bubble = brenth(to_solve, T_MAX, T_low_guess)
            except NotBoundedError:
                raise Exception("Bisecting solver could not find a solution between %g K and %g K" %(T_MAX, T_low_guess))

        return flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T_bubble, P=P)
    elif P is not None and VF is not None:
        T_low = flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, P=P, VF=1)[0]
        T_high = flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, P=P, VF=0)[0]
        info = []
        def err(T):
            T_calc, P_calc, VF_calc, xs, ys = flash_Tb_Tc_Pc(zs, Tbs, Tcs, Pcs, T=T, P=P)
            info[:] = T_calc, P_calc, VF_calc, xs, ys
            return VF_calc - VF
        P = brenth(err, T_low, T_high)
        return tuple(info)
    else:
        raise ValueError("Provide two of P, T, and VF")


def flash_ideal(zs, funcs, Tcs=None, T=None, P=None, VF=None):
    T_MAX = 50000.0
    N = len(zs)
    cmps = range(N)
    if Tcs is None: # numba: delete
        Tcs = [fi.solve_prop(1e6) for fi in funcs] # numba: delete
    if T is not None and P is not None:
        P_inv = 1.0/P
        Ks = [0.0]*N
        for i in cmps:
            Ks[i] = P_inv*funcs[i](T)
        ans = (T, P) + flash_inner_loop(zs=zs, Ks=Ks)
        return ans
    if T is not None and VF == 0.0:
        ys = [0.0]*N
        P_bubble = 0.0
        for i in cmps:
            v = funcs[i](T)*zs[i]
            P_bubble += v
            ys[i] = v

        P_inv = 1.0/P_bubble
        for i in cmps:
            ys[i] *= P_inv
        return (T, P_bubble, 0.0, zs, ys)
    if T is not None and VF == 1.0:
        xs = [0.0]*N
        P_dew = 0.
        for i in cmps:
            v = zs[i]/funcs[i](T)
            P_dew += v
            xs[i] = v
        P_dew = 1./P_dew
        for i in cmps:
            xs[i] *= P_dew
        return (T, P_dew, 1.0, xs, zs)
    elif T is not None and VF is not None:
        # Solve for in the middle of Pdew
        P_low = flash_ideal(zs, funcs, Tcs, T=T, VF=1)[1]
        P_high = flash_ideal(zs, funcs, Tcs, T=T, VF=0)[1]
        info = []
        def to_solve(P, info):
            T_calc, P_calc, VF_calc, xs, ys = flash_ideal(zs, funcs, Tcs, T=T, P=P)
            info[:] = T_calc, P_calc, VF_calc, xs, ys
            err = VF_calc - VF
            return err
        P = brenth(to_solve, P_low, P_high, args=(info,))
        return tuple(info)
    elif P is not None and VF == 1:
        def to_solve(T_guess):
            T_guess = abs(T_guess)
            P_dew = 0.
            for i in cmps:
                P_dew += zs[i]/funcs[i](T_guess)
            P_dew = 1./P_dew
            return P_dew - P

        # 2/3 average critical point
        T_guess = .66666*sum([Tcs[i]*zs[i] for i in cmps])
        try:
            T_dew = abs(secant(to_solve, T_guess, xtol=1e-12, maxiter=50))
        except Exception as e:
            T_dew = None
        if T_dew is None or T_dew > T_MAX*5.0:
            # Went insanely high T, bound it with brenth
            T_low_guess = sum([.1*Tcs[i]*zs[i] for i in cmps])
            bound = True
            try:
                err_low = to_solve(T_low_guess)
            except:
                bound = False
            try:
                err_high = to_solve(T_MAX)
            except:
                bound = False
            if bound and err_low*err_high > 0.0:
                bound = False

            if bound:
                T_dew = brenth(to_solve, T_low_guess, T_MAX, fa=err_low, fb=err_high)
            else:
                T_dew = secant(to_solve, min(min(Tcs)*0.9, T_guess), xtol=1e-12, maxiter=50, bisection=True, high=min(Tcs))
        
        xs = [P]*N
        for i in range(N):
            xs[i] *= zs[i]/funcs[i](T_dew)
        return (T_dew, P, 1.0, xs, zs)
    elif P is not None and VF == 0:
        def to_solve(T_guess):
            # T_guess = abs(T_guess)
            P_bubble = 0.0
            for i in cmps:
                P_bubble += zs[i]*funcs[i](T_guess)
            return P_bubble - P
        # 2/3 average critical point
        T_guess = sum([.55*Tcs[i]*zs[i] for i in cmps])
        try:
            T_bubble = abs(secant(to_solve, T_guess, maxiter=50, bisection=True, xtol=1e-12))
        except:
            T_bubble = None
        if T_bubble is None or T_bubble > T_MAX*5.0:
            # Went insanely high T, bound it with brenth
            T_low_guess = sum([.1*Tcs[i]*zs[i] for i in cmps])

            bound = True
            try:
                err_low = to_solve(T_low_guess)
            except:
                bound = False
            try:
                err_high = to_solve(T_MAX)
            except:
                bound = False
            if bound and err_low*err_high > 0.0:
                bound = False

            if bound:
                T_bubble = brenth(to_solve, T_low_guess, T_MAX, fa=err_low, fb=err_high)
            else:
                Tc_min = min(Tcs)
                T_bubble = secant(to_solve, min(Tc_min*0.9, T_guess), maxiter=50, bisection=True, high=Tc_min, xtol=1e-12)
        
        P_inv = 1.0/P
        ys = [0.0]*N
        for i in range(N):
            ys[i] = zs[i]*P_inv*funcs[i](T_bubble)
        return (T_bubble, P, 0.0, zs, ys)
    elif P is not None and VF is not None:
        bound = True
        try:
            T_low = flash_ideal(zs, funcs, Tcs, P=P, VF=1)[0]
            T_high = flash_ideal(zs, funcs, Tcs, P=P, VF=0)[0]
        except:
            bound = False
        info = []
        def err(T, zs, funcs, Tcs, P, VF, info, ignore_err):
            try:
                T_calc, P_calc, VF_calc, xs, ys = flash_ideal(zs, funcs, Tcs, T=T, P=P)
            except:
                if ignore_err:
                    return -0.5
                else:
                    raise ValueError("No solution in inner loop")
            info[:] = T_calc, P_calc, VF_calc, xs, ys
            return VF_calc - VF
        if bound:
            P = brenth(err, T_low, T_high, xtol=1e-14, args=(zs, funcs, Tcs, P, VF, info, False))
        else:
            T_guess = .5*sum([Tcs[i]*zs[i] for i in cmps])
            Tc_min = min(Tcs)
            # Starting at the lowest component's Tc should guarantee starting at two phases
            P = secant(err, Tc_min*(1.0-1e-7), xtol=1e-12, high=Tc_min, bisection=True,
                       args=(zs, funcs, Tcs, P, VF, info, True))
        return tuple(info)
    else:
        raise ValueError("Provide two of P, T, and VF")

