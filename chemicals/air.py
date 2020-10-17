# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020 Caleb Bell
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


.. contents:: :local:


"""
from __future__ import division
from fluids.numerics import trunc_exp
from math import exp, log, sqrt
__all__ = ['lemmon2000_air_A0', 'lemmon2000_air_dA0_dtau',
           'lemmon2000_air_d2A0_dtau2', 'lemmon2000_air_d3A0_dtau3',
           'lemmon2000_air_d4A0_dtau4',
           'lemmon2000_air_Ar', 'lemmon2000_air_dAr_dtau',
           'lemmon2000_air_d2Ar_dtau2', 'lemmon2000_air_d3Ar_dtau3',
           'lemmon2000_air_d4Ar_dtau4',
           'lemmon2000_air_dAr_ddelta', 'lemmon2000_air_d2Ar_ddelta2',
           'lemmon2000_air_d3Ar_ddelta3', 'lemmon2000_air_d4Ar_ddelta4']

# Get a good, fast variant of lemmon (2004) in here

# For values of tau above this, log(exp(87.31279*tau) + 2/3) reduces to 87.31279*tau in double precision
TAU_MAX_EXP_87 = 0.4207493606569795

def lemmon2000_air_A0(tau, delta):
    r'''Calculates the ideal gas Helmholtz energy of air according to Lemmon 
    (2000).
    
    .. math::
        \phi^\circ = \ln \delta + \sum_{i=1}^5 N_i\tau^{i-4} + N_6\tau^{1.5}
        + N_7\ln \tau + N_8\ln[1-\exp(-N_{11}\tau)] + N_9\ln[1-\exp(-N_{12}\tau)]
        + N_{10}\ln[2/3 + \exp(N_{13}\tau)]

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    A0 : float
        Ideal gas dimensionless Helmholtz energy A0/(RT) [-]

    Notes
    -----
            
    Examples
    --------
    >>> lemmon2000_air_A0(132.6312/200.0, 13000/10447.7)
    -14.65173785639
    '''
    tau_inv = 1.0/tau
    
#    exp0_00001 = exp(0.00001*tau)
    A0 =  (-0.00019536342*tau*sqrt(tau) + 17.275266575*tau + tau_inv*(tau_inv*(6.057194e-8*tau_inv 
            - 0.0000210274769) - 0.000158860716) + log(delta) + 2.490888032*log(tau)
    
            # These two logs both fail for tau < 1e-18, can be truncated but should not be necessary.
            + 0.791309509*log(1.0 - exp(-25.36365*tau)) + 0.212236768*log(1.0 - exp(-16.90741*tau)) 
            - 13.841928076)
    if tau < TAU_MAX_EXP_87:
        A0 -= 0.197938904*log(exp(87.31279*tau) + (2.0/3.0))
    else:
        A0 -= 17.282597957782162*tau # 17.282... = 87.31279*0.197938904
    return A0


def lemmon2000_air_dA0_dtau(tau, delta):
    r'''Calculates the first temperature derivative of ideal gas Helmholtz
    energy of air according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    dA0_dtau : float
        First derivative of `A0/RT` Ideal gas dimensionless Helmholtz energy
         with respect to `tau` [-]

    Notes
    -----
            
    Examples
    --------
    >>> lemmon2000_air_dA0_dtau(132.6312/200.0, 13000/10447.7)
    3.749095669249
    '''
    tau_inv = 1.0/tau
    dA0_dtau = (-0.00029304513*sqrt(tau) + tau_inv*(-tau_inv*(tau_inv*(1.8171582e-7*tau_inv 
            - 0.0000420549538) - 0.000158860716) + 2.490888032) + 17.275266575)
    try:
        # limit as tau gets high goes to zer0
        x0 = exp(87.31279*tau)
        dA0_dtau -= 17.282597957782162*x0/(x0 + (2.0/3.0))
    except:
        dA0_dtau -= 17.282597957782162
    try:
        # Limit is zero
        dA0_dtau += 20.0704974279478492/(exp(25.36365*tau) - 1.0)
    except:
        pass
    try:
        dA0_dtau += 3.58837405365087969/(exp(16.90741*tau) - 1.0)
    except:
        pass
    return dA0_dtau
    
    
def lemmon2000_air_d2A0_dtau2(tau, delta):
    r'''Calculates the second temperature derivative of ideal gas Helmholtz
    energy of air according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (126.192 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d2A0_dtau2 : float
        Second derivative of `A0/RT` Ideal gas dimensionless Helmholtz energy
         with respect to `tau` [-]

    Notes
    -----
            
    Examples
    --------
    >>> lemmon2000_air_d2A0_dtau2(132.6312/200.0, 13000/10447.7)
    -5.66675499015
    '''
    tau_inv = 1.0/tau
    d2A0_dtau2 = (-0.000146522565*sqrt(tau_inv) + tau_inv*(tau_inv*(tau_inv*(
            7.2686328e-7*tau_inv - 0.0001261648614) - 0.000317721432) 
            - 2.490888032)*tau_inv)
    
    if tau < 3.0:
        # 87.31279 Begins to have an impact a little under 0.5, others at 2.5 - set to 3 for safety
        x0 = exp(87.31279*tau)
        
        a = x0 + 2.0/3.0
        b = x0*(4.0/3.0 + x0) + 4.0/9.0
        d2A0_dtau2 += 1508.99184614226283*x0*(a*x0 - b)/(a*b)
                
        x0 = exp(25.36365*tau)
        a = (x0 - 1.0)
        b = x0*(x0 - 2.0) + 1.0
        d2A0_dtau2 -= 509.061072088369485*(a+b)/(a*b)
        
        x0 = exp(16.90741*tau)
        a = x0 - 1.0
        b = x0*(x0 - 2.0) + 1.0
        d2A0_dtau2 -= 60.670111358437417*(a+b)/(a*b)
    return d2A0_dtau2
    
def lemmon2000_air_d3A0_dtau3(tau, delta):
    r'''Calculates the third temperature derivative of ideal gas Helmholtz
    energy of air according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (126.192 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d3A0_dtau3 : float
        Third derivative of `A0/RT` Ideal gas dimensionless Helmholtz energy
         with respect to `tau` [-]

    Notes
    -----
            
    Examples
    --------
    >>> lemmon2000_air_d3A0_dtau3(132.6312/200.0, 13000/10447.7)
    17.10538866838
    '''
    tau_inv = 1.0/tau
    d3A0_dtau3 = (tau_inv*(0.0000732612825*sqrt(tau_inv) + tau_inv*(-tau_inv*(tau_inv*(
            3.6343164e-6*tau_inv - 0.0005046594456) - 0.000953164296) + 4.981776064)*tau_inv))
    
    if tau < 2.5:
        x0 = exp(16.90741*tau)
        x1 = exp(25.36365*tau)
        x3 = x0*x0
        
        x2 = exp(87.31279*tau)
        x4 = x2*x2
        x5 = x2*x2*x2
        d3A0_dtau3 += (-131754.288173931709*x2/(x2 + (2.0/3.0))
                + 395262.864521795127*x4/(1.333333333333333333*x2 + x4 + 0.444444444444444864) 
                - 263508.576347863418*x5/(1.333333333333333333*x2 + 2.0*x4 + x5 + 0.296296296296296668))
                
        d3A0_dtau3 += (25823.2937221483444/(3.0*x1 - 3.0*x1*x1 + x1*x1*x1 - 1.0) 
                            + 38734.9405832225166/(-2.0*x1 + x1*x1 + 1.0) 
                            + 12911.6468610741722/(x1 - 1.0))
        
        d3A0_dtau3 += (+ 2051.54889496551641/(3.0*x0 - 3.0*x3 + x0*x0*x0 - 1.0)
                        + 3077.32334244827462/(-2.0*x0 
                        + x3 + 1.0)  + 1025.77444748275821/(x0 - 1.0))
    return d3A0_dtau3

def lemmon2000_air_d4A0_dtau4(tau, delta):
    r'''Calculates the fourth temperature derivative of ideal gas Helmholtz
    energy of air according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d4A0_dtau4 : float
        Fourth derivative of `A0/RT` Ideal gas dimensionless Helmholtz energy
         with respect to `tau` [-]

    Notes
    -----
            
    Examples
    --------
    >>> lemmon2000_air_d4A0_dtau4(126.192/200.0, 13000/10447.7)
    -94.815532727
    '''
    tau_inv = 1.0/tau
    tau_inv2 = tau_inv*tau_inv
    d4A0_dtau4 = (-tau_inv2*(0.00010989192375*sqrt(tau_inv) - tau_inv2*(tau_inv*(
            tau_inv*(0.0000218058984*tau_inv - 0.002523297228)
            - 0.003812657184) - 14.945328192)))
    if tau < 0.4:
        x2 = exp(87.31279*tau)
        x5 = x2*x2
        x8 = x2*x2*x2
        x9 = x2*x2*x2*x2
        d4A0_dtau4 += 11503834.4949299842*x2*(-1.0/(x2 + 2.0/3.0) #1
        + x2*(7.0/((4.0/3.0)*x2 + x5 + (4.0/9.0))  #7
        - x2*(12.0/((4.0/3.0)*x2 + 2.0*x5 + x8 + (8.0/27.0)) #12
        + 6.0*x2/((-32.0/27.0)*x2 - (8.0/3.0)*(x5 +x8) - x9 - 16.0/81.0)  #6
        )))
    if tau < 2.0:
        x1 = exp(25.36365*tau)
        d4A0_dtau4 -= 327486.491907883901*(6.0/(x1*(x1*(x1*(x1 - 4.0) + 6.0) - 4.0) + 1.0) #6
        + 12.0/(x1*(x1*(x1 - 3.0) + 3.0) - 1.0) #12
        + 7.0/(x1*(x1 - 2.0) + 1.0) # 7
         + 1.0/(x1 - 1.0))
    if tau < 2.875:
        x0 = exp(16.90741*tau)
        d4A0_dtau4 += 17343.1891511144604*(6.0/(-x0*(x0*(x0*(x0 - 4.0) + 6.0) - 4.0) -1.0) #6
        - 12.0/(x0*(x0*(x0 - 3.0) + 3.0) - 1.0)  # 12
        - 7.0/(1.0 + x0*(x0 - 2.0)) # 7
        - 1.0/(x0 - 1.0)) # 1
    return d4A0_dtau4



def lemmon2000_air_Ar(tau, delta):
    r'''Calculates the residual Helmholtz energy of air according to Lemmon 
    (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    Ar : float
        Residual dimensionless Helmholtz energy Ar/(RT) [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt and many multiplies/adds.
            
    Examples
    --------
    >>> lemmon2000_air_Ar(132.6312/200.0, 13000/10447.7)
    -0.34683017661
    >>> lemmon2000_air_Ar(0.36842, 0.15880050154579475)
    0.0047988122806
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau12 = tau6*tau6
    tau_100 = tau**0.01
    tau2_100 = tau_100*tau_100
    tau4_100 = tau2_100*tau2_100
    tau5_100 = tau_100*tau4_100
    tau10_100 = tau5_100*tau5_100
    tau15_100 = tau5_100*tau10_100
    tau8_100 = tau4_100*tau4_100
    tau16_100 = tau8_100*tau8_100
    tau20_100 = tau4_100*tau16_100
    tau32_100 = tau16_100*tau16_100
    tau33_100 = tau_100*tau32_100
    tau64_100 = tau32_100*tau32_100
    tau80_100 = tau16_100*tau64_100
    tau40_100 = tau20_100*tau20_100
    tau97_100 = tau33_100*tau64_100
    tau45_100 = tau5_100*tau40_100
    tau90_100 = tau45_100*tau45_100
    tau160_100 = tau80_100*tau80_100
    tau320_100 = tau160_100*tau160_100
    x0 = exp(-delta)
    x1 = exp(-delta2)
    x2 = tau3*exp(-delta3)
    return (-0.101365037911999994*delta*tau160_100*x0 + 0.713116392079000017*delta*tau33_100
            - 0.146629609712999986*delta*tau40_100*x1*tau320_100 
            - 1.61824192067000006*delta*tau4_100*tau97_100 + 0.0148287891978000005*delta*taurt2*x2 
            + 0.118160747228999996*delta + 0.0714140178971000017*delta2 + 0.134211176704000013*delta3*tau15_100 
            - 0.031605587982100003*delta3*tau6*x1 - 0.17381369096999999*delta3*tau80_100*x0 
            - 0.00938782884667000057*delta3*x2*tau12 - 0.0865421396646000041*delta3 - 0.042053322884200002*delta4*tau20_100 
            + 0.0349008431981999989*delta4*tau2_100*tau33_100 + 0.0112626704218000001*delta4
            - 0.0472103183731000034*delta5*tau15_100*x0*tau80_100 + 0.000233594806141999996*delta5*tau3*taurt4*x1*delta6 
            - 0.0122523554252999996*delta6*tau*taurt4*x0 + 0.000164957183186000006*delta6*tau45_100*tau90_100)
    
    
def lemmon2000_air_dAr_dtau(tau, delta):
    r'''Calculates the first derivative of residual Helmholtz energy of air 
    with respect to tau according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    dAr_dtau : float
        First derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to tau, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt, 1 divisions
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_dAr_dtau(132.6312/200.0, 13000/10447.7)
    -1.8112257495223263
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau5 = tau*tau4
    tau10 = tau5*tau5
    tau_100 = tau**0.01
    tau2_100 = tau_100*tau_100
    tau4_100 = tau2_100*tau2_100
    tau8_100 = tau4_100*tau4_100
    tau16_100 = tau8_100*tau8_100
    tau32_100 = tau16_100*tau16_100
    tau33_100 = tau_100*tau32_100
    tau20_100 = tau4_100*tau16_100
    tau40_100 = tau20_100*tau20_100
    tau65_100 = tau32_100*tau33_100
    tau130_100 = tau65_100*tau65_100
    tau_inv_100 = 1.0/tau_100
    tau_inv2_100 = tau_inv_100*tau_inv_100
    tau_inv4_100 = tau_inv2_100*tau_inv2_100
    tau_inv5_100 = tau_inv_100*tau_inv4_100
    tau_inv8_100 = tau_inv4_100*tau_inv4_100
    tau_inv16_100 = tau_inv8_100*tau_inv8_100
    tau_inv20_100 = tau_inv4_100*tau_inv16_100
    tau_inv32_100 = tau_inv16_100*tau_inv16_100
    tau_inv64_100 = tau_inv32_100*tau_inv32_100
    tau_inv65_100 = tau_inv_100*tau_inv64_100
    tau_inv80_100 = tau_inv16_100*tau_inv64_100
    x0 = exp(-delta2)
    x1 = exp(-delta)
    x2 = delta6*taurt4
    x3 = exp(-delta3)
    return (-0.527866594966799996*delta*tau130_100*tau130_100*x0 + 0.0519007621922999984*delta*tau2*taurt2*x3 
            - 0.162184060659199991*delta*tau20_100*tau40_100*x1 - 1.63442433987669999*delta*tau_100 
            + 0.235328409386070025*delta*tau_inv2_100*tau_inv65_100 - 0.14081743270005001*delta3*tau10*tau4*x3 
            - 0.189633527892600018*delta3*tau5*x0 - 0.139050952775999992*delta3*tau_inv20_100*x1
            + 0.0201316765056000005*delta3*tau_inv5_100*tau_inv80_100 + 0.0122152951193699993*delta4*tau_inv65_100
            - 0.0084106645768400011*delta4*tau_inv80_100 + 0.000759183119961499907*delta5*tau2*x0*x2
            - 0.0448498024544450036*delta5*tau_inv5_100*x1 + 0.00022269219730110002*delta6*tau2_100*tau33_100 
            - 0.0153154442816249986*x1*x2)


def lemmon2000_air_d2Ar_dtau2(tau, delta):
    r'''Calculates the second derivative of residual Helmholtz energy of air 
    with respect to tau according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d2Ar_dtau2 : float
        Second derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to tau, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt, 2 divisions
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_d2Ar_dtau2(132.6312/200.0, 13000/10447.7)
    -0.7632109061747537
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    tau_inv = 1.0/tau
    taurt2 = sqrt(tau)
    tau_invrt2 = 1.0/taurt2
    tau_invrt4 = sqrt(tau_invrt2)
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau8 = tau4*tau4
    tau12 = tau4*tau8
    tau_inv_100 = tau_inv**0.01
    tau_inv2_100 = tau_inv_100*tau_inv_100
    tau_inv4_100 = tau_inv2_100*tau_inv2_100
    tau_inv8_100 = tau_inv4_100*tau_inv4_100
    tau_inv16_100 = tau_inv8_100*tau_inv8_100
    tau_inv32_100 = tau_inv16_100*tau_inv16_100
    tau_inv40_100 = tau_inv8_100*tau_inv32_100
    tau_inv33_100 = tau_inv_100*tau_inv32_100
    tau_inv65_100 = tau_inv32_100*tau_inv33_100
    tau_inv66_100 = tau_inv33_100*tau_inv33_100
    tau_inv99_100 = tau_inv33_100*tau_inv66_100
    tau_inv105_100 = tau_inv40_100*tau_inv65_100
    tau_inv80_100 = tau_inv40_100*tau_inv40_100
    tau_inv165_100 = tau_inv66_100*tau_inv99_100
    tau_inv20_100 = tau_inv4_100*tau_inv16_100
    tau_inv160_100 = tau_inv80_100*tau_inv80_100
    x0 = exp(-delta)
    x1 = tau_inv40_100*x0
    x2 = exp(-delta2)
    x3 = tau*exp(-delta3)
    return (-delta*(0.948167639463000089*delta2*tau4*x2 + 0.0171119250297600001*delta2*tau_inv80_100*tau_inv105_100
            - 0.0278101905551999921*delta2*x1*tau_inv80_100 + 1.9714440578007002*delta2*x3*tau12
            + 0.00793994182759050031*delta3*tau_inv165_100 - 0.00672853166147200157*delta3*tau_inv20_100*tau_inv160_100 
            - 0.00224249012272225226*delta4*tau_inv105_100*x0 - 0.00170816201991337503*delta5*tau*taurt2**0.5*x2*delta5
            - 0.0000779422690553850206*delta5*tau_inv65_100 + 0.00382886107040624965*delta5*tau_invrt2*tau_invrt4*x0 
            + 1.37245314691368003*tau**1.6*x2 + 0.15767003428866691*tau_inv2_100*tau_inv165_100 + 0.0163442433987670138*tau_inv99_100
            - 0.129751905480749996*taurt2*x3 + 0.0973104363955200058*x1))

def lemmon2000_air_d3Ar_dtau3(tau, delta):
    r'''Calculates the third derivative of residual Helmholtz energy of air 
    with respect to tau according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d3Ar_dtau3 : float
        Third derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to tau, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt, 4 divisions
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_d3Ar_dtau3(132.6312/200.0, 13000/10447.7)
    0.27922007457420
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    tau_inv = 1.0/tau
    taurt2 = sqrt(tau)
    tau_invrt2 = 1.0/taurt2
    tau_invrt4 = sqrt(tau_invrt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau_inv_100 = tau_inv**0.01
    tau_inv2_100 = tau_inv_100*tau_inv_100
    tau_inv4_100 = tau_inv2_100*tau_inv2_100
    tau_inv8_100 = tau_inv4_100*tau_inv4_100
    tau_inv16_100 = tau_inv8_100*tau_inv8_100
    tau_inv32_100 = tau_inv16_100*tau_inv16_100
    tau_inv33_100 = tau_inv_100*tau_inv32_100
    tau_inv66_100 = tau_inv33_100*tau_inv33_100
    tau_inv132_100 = tau_inv66_100*tau_inv66_100
    tau_inv140_100 = tau_inv8_100*tau_inv132_100
    tau_inv198_100 = tau_inv66_100*tau_inv132_100
    tau_inv41_100 = tau_inv8_100*tau_inv33_100
    tau_inv82_100 = tau_inv41_100*tau_inv41_100
    tau_inv164_100 = tau_inv82_100*tau_inv82_100
    tau_inv40_100 = tau_inv8_100*tau_inv32_100
    tau_inv44_100 = tau_inv4_100*tau_inv40_100
    tau_inv88_100 = tau_inv44_100*tau_inv44_100
    tau_inv264_100 = tau_inv132_100*tau_inv132_100
    tau_inv265_100 = tau_inv_100*tau_inv264_100
    tau_inv17_100 = tau_inv_100*tau_inv16_100
    tau_inv281_100 = tau_inv17_100*tau_inv264_100
    x0 = tau_inv132_100
    x1 = exp(-delta)
    x2 = exp(-delta3)
    x3 = exp(-delta2)
    return (delta*(-3.79267055785200036*delta2*tau3*x3 - 25.6287727514091017*delta2*tau6*x2*tau6 
           + 0.0316570613050560015*delta2*tau_inv4_100*tau_inv281_100 - 0.0333722286662399906*delta2*tau_inv88_100*x0*x1
           - 0.0121113569906496025*delta3*tau_inv140_100*tau_inv140_100 + 0.0131009040155243249*delta3*tau_inv265_100 
           - 0.00235461462885836479*delta4*tau_inv41_100*x1*tau_inv164_100
           + 0.00287164580280468724*delta5*tau_inv*tau_invrt2*tau_invrt4*x1 - 0.000050662474886000258*delta5*tau_inv33_100*x0 
           + 0.00213520252489171874*delta5*x3*delta5/tau_invrt4 - 2.19592503506188796*tau_inv2_100*tau_inv4_100/tau_inv66_100*x3 
           + 0.0389241745582079926*tau_inv140_100*x1 + 0.263308957262073706*tau_inv2_100*tau_inv265_100 
           + 0.0161808009647793419*tau_inv_100*tau_inv198_100 + 0.194627858221124994*taurt2*x2))
    
    
def lemmon2000_air_d4Ar_dtau4(tau, delta):
    r'''Calculates the fourth derivative of residual Helmholtz energy of air 
    with respect to tau according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d4Ar_dtau4 : float
        Fourth derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to tau, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt, 4 divisions
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_d4Ar_dtau4(132.6312/200.0, 13000/10447.7)
    -8.197368061417
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    tau_inv = 1.0/tau
    tau_invrt2 = sqrt(tau_inv)
    tau_invrt4 = sqrt(tau_invrt2)
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau8 = tau4*tau4
    tau10 = tau2*tau8
    tau_inv_100 = tau_inv**0.01
    tau_inv2_100 = tau_inv_100*tau_inv_100
    tau_inv4_100 = tau_inv2_100*tau_inv2_100
    tau_inv8_100 = tau_inv4_100*tau_inv4_100
    tau_inv16_100 = tau_inv8_100*tau_inv8_100
    tau_inv32_100 = tau_inv16_100*tau_inv16_100
    tau_inv40_100 = tau_inv8_100*tau_inv32_100
    tau_inv64_100 = tau_inv32_100*tau_inv32_100
    tau_inv80_100 = tau_inv16_100*tau_inv64_100
    tau_inv160_100 = tau_inv80_100*tau_inv80_100
    tau_inv240_100 = tau_inv80_100*tau_inv160_100
    tau_inv128_100 = tau_inv64_100*tau_inv64_100
    tau_inv256_100 = tau_inv128_100*tau_inv128_100
    tau_inv264_100 = tau_inv8_100*tau_inv256_100
    tau_inv265_100 = tau_inv_100*tau_inv264_100
    tau_inv34_100 = tau_inv2_100*tau_inv32_100
    tau_inv304_100 = tau_inv64_100*tau_inv240_100
    tau_inv320_100 = tau_inv64_100*tau_inv256_100
    tau_inv9_100 = tau_inv_100*tau_inv8_100
    tau_inv73_100 = tau_inv9_100*tau_inv64_100
    tau_inv146_100 = tau_inv73_100*tau_inv73_100
    tau_inv292_100 = tau_inv146_100*tau_inv146_100
    tau_inv68_100 = tau_inv34_100*tau_inv34_100
    tau_inv77_100 = tau_inv9_100*tau_inv68_100
    tau_inv145_100 = tau_inv68_100*tau_inv77_100
    tau_inv290_100 = tau_inv145_100*tau_inv145_100
    tau_inv20_100 = tau_inv4_100*tau_inv16_100
    tau_inv360_100 = tau_inv40_100*tau_inv320_100
    tau_inv384_100 = tau_inv128_100*tau_inv256_100
    x0 = exp(-delta)
    x1 = exp(-delta2)
    x2 = exp(-delta3)
    x3 = delta5*tau_invrt2*tau_invrt4
    return (-delta*(307.545273016909221*delta2*tau*x2*tau10 + 11.3780116735560011*delta2*tau2*x1 
            - 0.0734189030657279862*delta2*tau_inv320_100*x0 + 0.0902226247194096026*delta2*tau_inv_100*tau_inv384_100
            - 0.0339117995738188877*delta3*tau_inv20_100*tau_inv360_100 + 0.0347173956411394591*delta3*tau_inv73_100*tau_inv292_100 
            - 0.00482695998915964701*delta4*tau_inv_100*x0*tau_inv304_100 - 0.0000835930835619004249*delta5*tau_inv265_100
            + 0.00502538015490820288*tau_inv*x0*x3*tau_inv + 0.0544938443814911855*tau_inv240_100*x0 
            + 0.0321997939199108879*tau_inv34_100*tau_inv265_100 + 1.31755502103713296*tau_inv40_100*x1 
            + 0.703034915889736767*tau_inv77_100*tau_inv290_100 - 0.0973139291105624971*tau_invrt2*x2 
            - 0.000533800631222929685*x1*x3*delta5))

def lemmon2000_air_dAr_ddelta(tau, delta):
    r'''Calculates the first derivative of residual Helmholtz energy of air 
    with respect to delta according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    dAr_ddelta : float
        First derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to delta, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt,
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_dAr_ddelta(132.6312/200.0, 13000/10447.7)
    -1.8112257495223263
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau12 = tau6*tau6
    tau15 = tau3*tau12
    tau_100 = tau**0.01
    tau2_100 = tau_100*tau_100
    tau4_100 = tau2_100*tau2_100
    tau5_100 = tau_100*tau4_100
    tau10_100 = tau5_100*tau5_100
    tau15_100 = tau5_100*tau10_100
    tau8_100 = tau4_100*tau4_100
    tau16_100 = tau8_100*tau8_100
    tau20_100 = tau4_100*tau16_100
    tau32_100 = tau16_100*tau16_100
    tau33_100 = tau_100*tau32_100
    tau64_100 = tau32_100*tau32_100
    tau80_100 = tau16_100*tau64_100
    tau40_100 = tau20_100*tau20_100
    tau95_100 = tau15_100*tau80_100
    tau97_100 = tau33_100*tau64_100
    tau45_100 = tau5_100*tau40_100
    tau90_100 = tau45_100*tau45_100
    tau160_100 = tau80_100*tau80_100
    tau320_100 = tau160_100*tau160_100
    tau360_100 = tau40_100*tau320_100
    x0 = exp(-delta)
    x1 = 0.101365037911999994*tau160_100*x0
    x2 = exp(-delta2)
    x3 = tau360_100*x2
    x4 = exp(-delta3)
    x5 = 0.0281634865400100035*tau15*x4
    x6 = tau6*x2
    x7 = tau80_100*x0
    x8 = tau95_100*x0
    x9 = tau3*taurt2*x4
    x10 = tau*x0
    x11 = delta6*taurt4
    x12 = tau3*x2
    x13 = delta6
    return (delta*x1 + 0.142828035794200003*delta + 0.402633530112000038*delta2*tau15_100 
            + 0.293259219425999973*delta2*x3 - delta2*x5 - 0.0948167639463000089*delta2*x6 
            - 0.52144107290999997*delta2*x7 - 0.259626418993800012*delta2 - 0.168213291536800008*delta3*tau20_100
            + 0.139603372792799996*delta3*tau2_100*tau33_100 + 0.17381369096999999*delta3*x7 
            - 0.0444863675934000016*delta3*x9 + 0.0450506816872000004*delta3
            + 0.00256954286756199985*delta4*taurt4*x12*x13 + 0.063211175964200006*delta4*x6
            - 0.23605159186550001*delta4*x8 + 0.000989743099116000089*delta5*tau45_100*tau90_100 
            - 0.0735141325518000044*delta5*taurt4*x10 + delta5*x5 + 0.0472103183731000034*delta5*x8 
            + 0.713116392079000017*tau33_100 - 1.61824192067000006*tau4_100*tau97_100 - x1
            + 0.0122523554252999996*x10*x11 - 0.000467189612283999993*x11*x12*x13 - 0.146629609712999986*x3 
            + 0.0148287891978000005*x9 + 0.118160747228999996)
    
    
def lemmon2000_air_d2Ar_ddelta2(tau, delta):
    r'''Calculates the second derivative of residual Helmholtz energy of air 
    with respect to delta according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d2Ar_ddelta2 : float
        Second derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to delta, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt,
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_d2Ar_ddelta2(132.6312/200.0, 13000/10447.7)
    0.27027259528316
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    delta7 = delta*delta6
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau12 = tau6*tau6
    tau15 = tau3*tau12
    tau_20 = tau**0.05
    tau2_20 = tau_20*tau_20
    tau3_20 = tau_20*tau2_20
    tau4_20 = tau2_20*tau2_20
    tau8_20 = tau4_20*tau4_20
    tau16_20 = tau8_20*tau8_20
    tau18_20 = tau2_20*tau16_20
    tau19_20 = tau_20*tau18_20
    tau9_20 = tau_20*tau8_20
    tau32_20 = tau16_20*tau16_20
    tau64_20 = tau32_20*tau32_20
    tau72_20 = tau8_20*tau64_20
    x0 = exp(-delta)
    x1 = tau32_20*x0
    x2 = exp(-delta3)
    x3 = tau15*x2
    x4 = 1.04288214581999994*tau16_20*x0
    x5 = exp(-delta2)
    x6 = tau6*x5
    x7 = tau72_20*x5
    x8 = delta3*x0
    x9 = tau19_20*x0
    x10 = tau3*taurt2*x2
    x11 = tau*x0
    x12 = taurt4*x11
    x13 = delta6*taurt4
    x14 = tau3*x5
    x15 = delta4*taurt4*x14
    x16 = delta7
    return (0.805267060224000075*delta*tau3_20 - 0.101365037911999994*delta*x1 - 0.0563269730800200069*delta*x3 
            - delta*x4 - 0.189633527892600018*delta*x6 + 0.879777658277999919*delta*x7 
            - 0.519252837987600024*delta + 0.418810118378399987*delta2*tau3_20*tau4_20 
            - 0.504639874610400052*delta2*tau4_20 - 0.177945470373600007*delta2*x10 + delta2*x4 
            + 0.135152045061600001*delta2 + 0.442478231749400042*delta3*x6 - 0.586518438851999946*delta3*x7 
            + 0.00494871549558000001*delta4*tau9_20*tau18_20 - 0.367570662759000022*delta4*x12
            + 0.225307892320080028*delta4*x3 + 0.47210318373100002*delta4*x9 + 0.133459102780200012*delta5*x10 
            + 0.147028265103600009*delta5*x12 - 0.126422351928400012*delta5*x6 - 0.0472103183731000034*delta5*x9 
            - 0.0844904596200300173*delta7*x3 - 0.17381369096999999*tau16_20*x8 - 0.94420636746200004*tau19_20*x8 
            + 0.202730075823999989*x1 - 0.0122523554252999996*x11*x13 + 0.000934379224567999985*x13*x14*x16 
            - 0.0107453610825320005*x15*x16 + 0.0256954286756199968*x15*delta5 + 0.142828035794200003)

def lemmon2000_air_d3Ar_ddelta3(tau, delta):
    r'''Calculates the third derivative of residual Helmholtz energy of air 
    with respect to delta according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d3Ar_ddelta3 : float
        Third derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to delta, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt,
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_d3Ar_ddelta3(132.6312/200.0, 13000/10447.7)
    0.1849386546766
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    delta8 = delta4*delta4
    delta12 = delta4*delta8
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau12 = tau6*tau6
    tau15 = tau3*tau12
    tau_20 = tau**0.05
    tau2_20 = tau_20*tau_20
    tau3_20 = tau_20*tau2_20
    tau4_20 = tau2_20*tau2_20
    tau8_20 = tau4_20*tau4_20
    tau16_20 = tau8_20*tau8_20
    tau18_20 = tau2_20*tau16_20
    tau19_20 = tau_20*tau18_20
    tau9_20 = tau_20*tau8_20
    tau32_20 = tau16_20*tau16_20
    tau64_20 = tau32_20*tau32_20
    tau72_20 = tau8_20*tau64_20
    x0 = exp(-delta3)
    x1 = tau15*x0
    x2 = exp(-delta)
    x3 = tau16_20*x2
    x4 = tau32_20*x2
    x5 = exp(-delta2)
    x6 = tau6*x5
    x7 = tau72_20*x5
    x8 = tau19_20*x2
    x9 = 2.83261910238600034*x8
    x10 = tau3*taurt2*x0
    x11 = delta*x10
    x12 = tau3*taurt4*x5
    x13 = tau*taurt4*x2
    x14 = delta8
    x15 = delta2*x12
    return (0.837620236756799974*delta*tau3_20*tau4_20 - 1.0092797492208001*delta*tau4_20 
            + 0.25347137886009008*delta*x1*x14 + 3.1286464374599996*delta*x3 + 0.101365037911999994*delta*x4
            + 0.270304090123200003*delta + 0.0336376520844480012*delta12*x12 - 1.5643232187299998*delta2*x3 
            + 1.70670175103340016*delta2*x6 - 3.51911063311199968*delta2*x7 - delta2*x9
            + 0.01979486198232*delta3*tau9_20*tau18_20 + 1.07021248852038009*delta3*x1 
            - 1.47028265103600009*delta3*x13 + 0.17381369096999999*delta3*x3 + delta3*x9 
            + 1.20113192502180022*delta4*x10 + 1.10271198827700001*delta4*x13 - 1.51706822314080014*delta4*x6
            + 1.17303687770399989*delta4*x7 - 0.708154775596500086*delta4*x8 - 0.220542397655400013*delta5*x13 
            + 0.0472103183731000034*delta5*x8 - 1.26735689430045029*delta6*x1 + 0.0122523554252999996*delta6*x13 
            + 0.252844703856800024*delta6*x6 + 0.231258858080579971*delta8*x12 + 0.805267060224000075*tau3_20 
            - 0.0563269730800200069*x1 - 0.355890940747200013*x11 - 0.400377308340600035*x11*delta6 
            - 0.169589829259091995*x14*x15 - 0.00186875844913599997*x15*delta12 - 1.04288214581999994*x3 
            - 0.304095113735999956*x4 - 0.189633527892600018*x6 + 0.879777658277999919*x7 - 0.519252837987600024)

def lemmon2000_air_d4Ar_ddelta4(tau, delta):
    r'''Calculates the fourth derivative of residual Helmholtz energy of air 
    with respect to delta according to Lemmon (2000).
    
    Parameters
    ----------
    tau : float
        Dimensionless temperature, (132.6312 K)/T [-]
    delta : float
        Dimensionless density, rho/(10447.7 mol/m^3), [-]

    Returns
    -------
    d4Ar_ddelta4 : float
        Fourth derivative of residual dimensionless Helmholtz energy Ar/(RT)
        with respect to delta, [-]

    Notes
    -----
    The cost of this function is 1 power, 3 exp, 2 sqrt,
    and the necessary adds/multiplies.
            
    Examples
    --------
    >>> lemmon2000_air_d4Ar_ddelta4(132.6312/200.0, 13000/10447.7)
    0.37902213262258
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta2*delta4
    delta7 = delta*delta6
    delta8 = delta4*delta4
    delta9 = delta*delta8
    delta11 = delta2*delta9
    taurt2 = sqrt(tau)
    taurt4 = sqrt(taurt2)
    tau2 = tau*tau
    tau3 = tau*tau2
    tau6 = tau3*tau3
    tau12 = tau6*tau6
    tau15 = tau3*tau12
    tau_20 = tau**0.05
    tau2_20 = tau_20*tau_20
    tau4_20 = tau2_20*tau2_20
    tau6_20 = tau2_20*tau4_20
    tau8_20 = tau4_20*tau4_20
    tau16_20 = tau8_20*tau8_20
    tau18_20 = tau2_20*tau16_20
    tau19_20 = tau_20*tau18_20
    tau9_20 = tau_20*tau8_20
    tau32_20 = tau16_20*tau16_20
    tau64_20 = tau32_20*tau32_20
    tau72_20 = tau8_20*tau64_20
    x0 = exp(-delta)
    x1 = tau16_20*x0
    x2 = tau32_20*x0
    x3 = tau19_20*x0
    x4 = 5.66523820477200069*x3
    x5 = exp(-delta2)
    x6 = tau6*x5
    x7 = tau72_20*x5
    x8 = exp(-delta3)
    x9 = tau15*x8
    x10 = tau3*taurt2*x8
    x11 = tau3*x5
    x12 = taurt4*x11
    x13 = tau*x0
    x14 = taurt4*x13
    x15 = delta4*taurt4
    x16 = delta9
    return (-6.25729287491999919*delta*x1 - 0.101365037911999994*delta*x2 - delta*x4 
            + 3.79267055785200036*delta*x6 - 8.79777658277999919*delta*x7 + 0.742831483531560033*delta11*x12
            - 0.760414136580270239*delta11*x9 + 0.0593845859469600001*delta2*tau9_20*tau18_20 
            + 2.08576429163999988*delta2*x1 - 4.41084795310800004*delta2*x14 + 11.3304764095440014*delta2*x3
            + 3.37961838480120003*delta2*x9 - 0.17381369096999999*delta3*x1 + 5.87220052232880096*delta3*x10
            + 5.88113060414400035*delta3*x14 - delta3*x4 - 9.48167639463*delta3*x6 
            + 11.7303687770399989*delta3*x7 + 0.94420636746200004*delta4*x3 + 0.294056530207200018*delta5*x14 
            - 0.0472103183731000034*delta5*x3 + 4.55120466942240043*delta5*x6 - 2.34607375540799978*delta5*x7
            - 10.8147788313638422*delta5*x9 - 6.40603693344960057*delta6*x10 + 0.00373751689827199994*delta6*x12*x16 
            - 0.0122523554252999996*delta6*x14 + 1.85007086464463977*delta7*x12
            - 0.505689407713600048*delta7*x6 + 6.08331309264216191*delta8*x9 + 1.20113192502180022*delta9*x10 
            - 2.15841600875207984*delta9*x12 - 1.0092797492208001*tau4_20 + 0.837620236756799974*tau_20*tau6_20 
            + 4.17152858327999976*x1 - 0.355890940747200013*x10 - 0.0934379224567999933*x11*x15*x16 
            - 2.20542397655400002*x13*x15 + 0.405460151647999978*x2 + 0.270304090123200003)
