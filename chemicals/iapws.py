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

"""
from __future__ import division
from math import exp, log, sqrt

__all__ = ['iapws97_dG_dpi_region1', 'iapws97_dGr_dpi_region2', 'iapws97_dGr_dpi_region5',
           
           'iapws97_boundary_2_3', 'iapws97_boundary_3uv', 'iapws97_boundary_3ef', 
           'iapws97_boundary_3ef', 'iapws97_boundary_3cd', 'iapws97_boundary_3gh',
           'iapws97_boundary_3ij', 'iapws97_boundary_3jk', 'iapws97_boundary_3mn', 
           'iapws97_boundary_3qu', 'iapws97_boundary_3rx', 'iapws97_boundary_3wx',
           'iapws97_boundary_3ab', 'iapws97_boundary_3op'
           
           ]


# The intention is to have a partial sovler for IAPWS-95 and IAPWS-97
# Which does not compute properties in general; the standard has dozens
# and dozens more are needed by a thermodynamic package

# I am just looking to get as fast as possible solvers for:
# IAPWS-97: T, P -> rho
# IAPWS-95: T, P -> rho using IAPWS-97 as initial guess
# IAPWS-95: T, rho -> P
# IAPWS-95: P, rho -> V
# I have been casually working this for over 5 years


def iapws97_boundary_2_3(T):
    '''Above this pressure we are in region 3.
    
    >>> iapws97_boundary_2_3(0.623150000E3)
    16529164.2526216
    '''
    return (0.34805185628969E9 - T*(0.11671859879975E7 - 0.10192970039326E4*T))

def iapws97_boundary_3uv(P):
    '''
    >>> iapws97_boundary_3uv(22.3E6)
    647.7996121480069'''
    return (P*(P*(2.867916822636969863e-21*P - 2.228141349037550121e-13) 
               + 8.905796021353068107e-6) + 528.1996462630620499)
#    P = P/1E6
#    T = sum([nis3uv[i]*P**Iis3uv[i] for i in range(4)])
#    return T

def iapws97_boundary_3ef(P):
    '''
    >>> iapws97_boundary_3ef(40E6)
    713.959399239744
    '''
    return 3.7278880039999996e-6*P + 564.843879079744056
#    P = P/1E6
#    T = 3.727888004*(P-22.064)+647.096
#    return T

def iapws97_boundary_3cd(P):
    '''
    >>> iapws97_boundary_3cd(25E6)
    649.3659208321279'''
    return (P*(P*(1.59090746562728991e-22*P - 1.2728354929587799e-14) 
               + 2.78233532206914969e-6) + 585.27696669634895)
#    P = P/1E6
#    T = sum([nis3cd[i]*P**Iis3cd[i] for i in range(4)])
#    return T

def iapws97_boundary_3gh(P):
    '''
    >>> iapws97_boundary_3gh(25E6)
    656.69805722612
    '''
    return (P*(P*(P*(7.51608051114156857e-18 - 7.87105249910382769e-26*P)
                  - 2.69029173140129987e-10) + 0.0042814358479154593) - 24928.4240900418008)
#    P = P/1E6
#    T = sum([nis3gh[i]*P**Iis3gh[i] for i in range(5)])
#    return T

def iapws97_boundary_3ij(P):
    '''
    >>> iapws97_boundary_3ij(25E6)
    660.7865756716819'''
    return (P*(P*(P*(5.15308185433081825e-29*P - 5.87071076864458977e-21) 
                  + 2.60763050899561981e-13) - 6.16179320924617007e-7) + 584.814781649163024)
#    P = P/1E6
#    T = sum([nis3ij[i]*P**Iis3ij[i] for i in range(5)])
#    return T

def iapws97_boundary_3jk(P):
    '''
    >>> iapws97_boundary_3jk(25E6)
    668.1915358826951'''
    return (P*(P*(P*(1.37897492684193974e-28*P - 1.57391839848015003e-20) 
                  + 6.97072596851896056e-13) - 7.70600270141674947e-6) + 617.229772068439047)
#    P = P/1E6
#    T = sum([nis3jk[i]*P**Iis3jk[i] for i in range(5)])
#    return T

def iapws97_boundary_3mn(P):
    '''
    >>> iapws97_boundary_3mn(22.8E6)
    649.6054132953997'''
    return (P*(P*(1.92871054508107992e-21*P - 1.58365725441647998e-13)
               + 7.61978122720127966e-6) + 535.339483742384004)
#    P = P/1E6
#    T = sum([nis3mn[i]*P**Iis3mn[i] for i in range(4)])
#    return T

def iapws97_boundary_3qu(P):
    '''
    >>> iapws97_boundary_3qu(22E6)
    645.6355027340121'''
    return (P*(P*(1.22240301070144985e-21*P - 1.02020639611015996e-13) 
               + 5.29062258221221963e-6) + 565.60364823912596)
#    P = P/1E6
#    T = sum([nis3qu[i]*P**Iis3qu[i] for i in range(4)])
#    return T

def iapws97_boundary_3rx(P):
    '''
    >>> iapws97_boundary_3rx(22E6)
    648.26227536701
    '''
    return (P*(P*(2.43293362700451993e-13 - 2.94905044740798962e-21*P)
               - 1.02961025163668997e-6) + 584.561202520005963)
#    P = P/1E6
#    T = sum([nis3rx[i]*P**Iis3rx[i] for i in range(4)])
#    return T

### Log P boundaries
def iapws97_boundary_3wx(logP_MPa, logP_MPa_inv):
    '''
    >>> iapws97_boundary_3wx(log(22.3), 1/log(22.3))
    648.204947950734
    '''
    return (logP_MPa*(14.7370491183190993*logP_MPa + 97.3505869861951965) 
            + 7.28052609145380014 + logP_MPa_inv*(329.196213998375015 + 873.371668682416953*logP_MPa_inv))
#    P = P/1E6
#    T = sum([nis3wx[i]*logP_MPa**Iis3wx[i] for i in range(5)])
#    return T

def iapws97_boundary_3ab(logP_MPa, logP_MPa_inv):
    '''
    >>> iapws97_boundary_3ab(log(40), 1/log(40))
    693.0341408296053'''
    return (logP_MPa*(21.3144632222113017*logP_MPa - 187.661219490113012)
            + 1547.93642129415002 - logP_MPa_inv*(1918.87498864292002 - 918.419702359447001*logP_MPa_inv))
#    T = sum([nis3ab[i]*logP_MPa**Iis3ab[i] for i in range(5)])
#    return T

def iapws97_boundary_3op(logP_MPa, logP_MPa_inv):
    '''
    >>> iapws97_boundary_3op(log(22.8), 1/log(22.8))
    650.010694314133'''
    return (logP_MPa*(64.2859598466067013*logP_MPa - 332.500170441277987) 
            + 969.461372400213008 + logP_MPa_inv*(773.845935768222034 - 1523.13732937084001*logP_MPa_inv))
#    T = sum([nis3op[i]*logP_MPa**Iis3op[i] for i in range(5)])
#    return T




def iapws97_dG_dpi_region1(tau, pi):
    r'''Calculates dG_dpi for region 1.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (1386 K)/T [-]
    pi : float
        Dimensionless pressure, P/(16.53 MPa), [-]

    Returns
    -------
    dG_dpi : float
        Derivative of dimentionless Gibbs energy G/(RT) with respect to `pi`,
        [-]

    Notes
    -----
    Used in density solution.

    Examples
    --------
    >>> iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6)
    0.1292327182544
    '''
    pit = 7.1 - pi
    taut = tau - 1.222
    pit2 = pit*pit
    pit3 = pit2*pit
    pit13 = pit3*pit3
    pit13 *= pit13*pit
    taut2 = taut*taut
    taut3 = taut2*taut
    taut4 = taut2*taut2
    taut8 = taut4*taut4
    taut_inv = 1.0/taut
    taut_inv2 = taut_inv*taut_inv
    taut_inv3 = taut_inv2*taut_inv
    taut_inv4 = taut_inv2*taut_inv2
    taut_inv8 = taut_inv4*taut_inv4
    taut_inv32 = taut_inv8*taut_inv8
    taut_inv32 *= taut_inv32
    
    return (pit*(pit*(pit*(pit*(pit3*(pit13*(pit2*(pit3*pit3*(pit*(pit*(
            taut_inv32*taut_inv8*(2.99318679335865594e-24*pit*taut_inv - 5.65070932023524029e-23))
    + 3.58428679202129995e-22*taut_inv32*taut_inv8*taut) 
    - 7.6373766822105502e-22*taut_inv32*taut_inv3*taut_inv3) - 3.33001080055983015e-19*taut_inv32*taut)
    + 1.44400475720615078e-17*taut_inv32*taut_inv3) + 1.0*(1.39398969845072005e-9*taut4*taut
    + 1.01874413933127995e-8)*taut_inv8*taut_inv3) + 2.02584984300584983e-6*taut_inv8) 
    + 1.0*(taut3*(5.73669197516959969e-13*taut8*taut4 + 2.60684891582404009e-6)
    + 8.97011276319999943e-6)*taut_inv4*taut_inv) + (taut4*(2.55615384360309023e-9*taut4*taut2 
    + 8.48123939559359992e-6) + 0.0000950389345351620047)*taut_inv4)
    + (taut3*(1.45389992595188003e-15*taut8*taut8*taut + 8.82836906616919934e-6*taut3
   - 0.0000953227878139740028*taut + 0.000600035615860519981)
   + 0.000943686421465339954)*taut_inv3) + (0.0000528383579699300023*taut8*taut4
   + taut8*((0.0218417171754139994*taut + 0.0325297487705049973)*taut + 0.0189900682184190005)
   + 0.000607063015658739955*taut2 - 0.000283190801238040004)*taut_inv8*taut_inv)


def iapws97_dGr_dpi_region2(tau, pi):
    r'''Calculates dGr_dpi for region 2.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (540 K)/T [-]
    pi : float
        Dimensionless pressure, P/(1 MPa), [-]

    Returns
    -------
    dGr_dpi : float
        Derivative of dimentionless residual Gibbs energy G/(RT) with respect 
        to `pi`, [-]

    Notes
    -----
    Used in density solution.

    Examples
    --------
    >>> iapws97_dGr_dpi_region2(.656, 16)
    -0.0062926319312
    '''
    taut = tau - 0.5
    pi2 = pi*pi
    taut2 = taut*taut
    taut3 = taut*taut2
    taut4 = taut2*taut2
    taut6 = taut4*taut2
    taut8 = taut4*taut4
    taut13 = taut6*taut4*taut3
    taut21 = taut13*taut8
    taut29 = taut21*taut8
    # 53 from 13*13*!3*!3*3
    # 57 from 
    return (pi*(pi*(pi*(pi*(pi*(pi*(pi*(pi*(pi*(pi2*pi2*pi2*(pi2*(pi2
        *(pi*(pi*(pi*(pi*taut13*taut13*(taut13*taut*(1.32995316841867198e-15 - 0.0000226487297378903988*taut13*taut4*taut)
        + 1.75410265428146418e-27) - 2.93678005497663003e-14*taut29*taut8*taut2) + 0.0000832192847496054092*taut*taut21*taut29*taut3)
    - 1.24017662339841913e-24*taut21) + taut8*taut8*taut4*(taut13*taut2*(6.1258633752463995e-12
    - 0.0000840049353964159951*taut13) + 1.78371690710842005e-23)) - 6.05920510335077989*taut21*taut21*taut13*taut2) 
    + taut29*(1.71088510070543998*taut21 - 1.29412653835175996e-9)) + taut4*(taut6
    *(-1.00181793795109993e-8*taut4 - 1.02347470959289996e-12) + 1.04069652101739995e-18))
    + 1.78287415218792009e-7*taut13) + taut4*taut4*(9.00496908836719986e-11 - 65.8490727183984035*taut13*taut13*taut2))
    + taut6*taut4*taut*(-0.27262789705017304*taut6*taut6*taut2 - 8.83526622937069987e-6) - 4.13416950269890026e-17)
    + taut3*(taut13*(-143.374451604623999*taut13*taut6 - 0.012702883392812999) - 1.00288598706366e-10))
    + 0.0000114610381688305001*taut6*taut) + taut*(taut*(1.92901490874028006e-6*taut + 5.11628714091400033e-8)
    - 3.15389238237468004e-9)) + taut*(taut2*(taut3*(-0.122004760687946995*taut21*taut8 
    - 0.00451017736264439952) - 0.0000968330317157100001) + 1.31612001853305008e-6) + 6.14452130769269999e-8)
    + taut*(taut*(taut2*(taut3*(-0.0000533490958281740028*taut21*taut8 - 0.0875945913011459965) 
    - 0.00787855544867100029) - 0.000378979750326299998) - 0.000066065283340406)) 
    + taut*(taut*(taut*(-0.0503252787279300021*taut3 - 0.0575812590834320001) - 0.0459960136963650026) 
    - 0.0178348622923579989) - 0.00177317424732129992)


def iapws97_dGr_dpi_region5(tau, pi):
    r'''Calculates dGr_dpi for region 5.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (1000 K)/T [-]
    pi : float
        Dimensionless pressure, P/(1 MPa), [-]

    Returns
    -------
    dGr_dpi : float
        Derivative of dimentionless residual Gibbs energy G/(RT) with respect 
        to `pi`, [-]

    Notes
    -----
    Used in density solution.

    Examples
    --------
    >>> iapws97_dGr_dpi_region5(.5, 30.0)
    0.0004009761854002
    '''
    tau3 = tau*tau
    tau3 *= tau
    return (pi*tau3*(4.48800748189699983e-6 + tau3*(1.13758364468865005e-7*pi*tau - 8.23265509069419973e-6*tau3))
            + tau*(tau*(0.000901537616739440007 - 0.00502700776776479966*tau) + 0.00157364048552589993))
