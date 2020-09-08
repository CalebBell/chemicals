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
from chemicals.vapor_pressure import Psat_IAPWS, Tsat_IAPWS

__all__ = ['iapws97_dGr_dpi_region2', 'iapws97_dGr_dpi_region5',
           'iapws97_boundary_2_3', 'iapws97_boundary_3uv', 'iapws97_boundary_3ef', 
           'iapws97_boundary_3ef', 'iapws97_boundary_3cd', 'iapws97_boundary_3gh',
           'iapws97_boundary_3ij', 'iapws97_boundary_3jk', 'iapws97_boundary_3mn', 
           'iapws97_boundary_3qu', 'iapws97_boundary_3rx', 'iapws97_boundary_3wx',
           'iapws97_boundary_3ab', 'iapws97_boundary_3op',
           'iapws97_identify_region_TP', 'iapws97_region_3', 'iapws97_region3_rho',
           'iapws97_rho',
           
           'iapws95_d2A_d2deltar', 'iapws95_dA_ddeltar',
           
           
           'iapws97_G_region1', 'iapws97_dG_dpi_region1', 'iapws97_d2G_dpi2_region1',
           ]

__numba_additional_funcs__ = ['iapws97_region3_a', 'iapws97_region3_b', 'iapws97_region3_c', 
    'iapws97_region3_d', 'iapws97_region3_e', 'iapws97_region3_f', 'iapws97_region3_g', 'iapws97_region3_h',
    'iapws97_region3_i', 'iapws97_region3_j', 'iapws97_region3_k', 'iapws97_region3_l', 'iapws97_region3_m',
    'iapws97_region3_n', 'iapws97_region3_o', 'iapws97_region3_p', 'iapws97_region3_q', 'iapws97_region3_r', 
    'iapws97_region3_s', 'iapws97_region3_t', 'iapws97_region3_u', 'iapws97_region3_v', 'iapws97_region3_w', 
    'iapws97_region3_x', 'iapws97_region3_y', 'iapws97_region3_z']



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


### Fast dG_dpi and dGr_dpi for density calls
def iapws97_G_region1(tau, pi):
    pit = 7.1 - pi
    taut = tau - 1.222
    taut_inv = 1.0/taut
    
    pit2 = pit*pit
    pit4 = pit2*pit2
    taut_inv2 = pit4*pit4 # abuse taut_inv2 variable as a temporary
    pit21 = taut_inv2*taut_inv2*pit4*pit
    pit29 = pit21*taut_inv2
    
    taut_inv2 = taut*taut # abuse taut_inv2 variable as a temporary
    taut3 = taut_inv2*taut
    taut5 = taut_inv2*taut3
    
    taut_inv2 = taut_inv*taut_inv
    taut_inv3 = taut_inv2*taut_inv
    taut_inv9 = taut_inv3*taut_inv3*taut_inv3
    taut_inv29 = taut_inv9*taut_inv9*taut_inv9*taut_inv2
    return (pit*(-0.02184171717541399937*taut - 0.01899006821841900047*taut_inv 
            - 0.0325297487705049973 - 0.00005283835796993000233*taut3 
            - 0.000607063015658739955*taut_inv3*taut_inv3*taut_inv 
            + 0.0002831908012380400042*taut_inv9) 
        
           + taut*(3.385516916838500201
            + 0.00004766139390698700138*pit2 
            - 0.9579196338787200338*taut)
            - 0.8454818716911399745*taut_inv
            - 3.756360367204000017
            + 0.1463297121316700089*taut_inv2 
            
            + pit2*(-0.0003000178079302599906
            - 4.414184533084599669e-6*taut3 - 0.0004718432107326699771*taut_inv3 
            - 7.269499629759400146e-16*taut5*taut5*taut5*taut*taut 
            + pit*(-2.827079798531199973e-6 - 0.00003167964484505400157*taut_inv3*taut_inv 
            - 8.520512812010300437e-10*taut5*taut))
            
            + taut3*(0.157720385132280011 - 0.01661641719950100043*taut)
            + 0.0008121462998356799657*taut5
            
            + pit4*(taut_inv2*(- 6.517122289560100218e-7 - taut_inv3*(2.242528190799999857e-6 
            + 4.051699686011699983e-7*pit*taut_inv3)
            - 1.273430174164099942e-9*pit4*taut_inv9)
            - 1.434172993792399922e-13*taut5*taut5 
            - 1.742487123063400057e-10*pit4*taut_inv3*taut_inv3)

            + taut_inv29*(1.44783078285210013e-20*pit21*pit2*taut_inv2 
            - 6.876213129553099646e-19*pit21 
            + pit29*taut_inv9*(2.633578166279499979e-23
            - 1.194762264007099993e-23*pit*taut_inv 
            + 1.822809458140400033e-24*pit2*taut_inv2
            - 9.353708729245799802e-26*pit2*pit*taut_inv3)))

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
    pit7 = pit3*pit3*pit
    pit14 = pit7*pit7
    taut3 = taut*taut*taut
    taut6 = taut3*taut3
    taut10 = taut6*taut3*taut
    taut_inv = 1.0/taut
    taut_inv2 = taut_inv*taut_inv
    taut_inv4 = taut_inv2*taut_inv2
    taut_inv5 = taut_inv4*taut_inv
    taut_inv7 = taut_inv5*taut_inv2
    taut_inv17 = taut_inv5*taut_inv5*taut_inv7
    pit14taut_inv17 = taut_inv17*pit14
    return (pit*(-0.0000953227878139740028*taut + 0.000600035615860519981 + 8.82836906616919934e-6*(taut3) 
            + 0.000943686421465339954*(taut_inv2*taut_inv) + 1.45389992595188003e-15*(taut10*taut6*taut))
            + 0.0218417171754139994*taut + 0.0189900682184190005*taut_inv + 0.0325297487705049973 
            + pit2*(8.48123939559359992e-6 + 0.0000950389345351620047*(taut_inv4) +
            2.55615384360309023e-9*(taut6)) + 2.60684891582404009e-6*(pit3)*(taut_inv2)
            + 8.97011276319999943e-6*(pit3)*(taut_inv5) + 5.73669197516959969e-13*(pit3)*(taut10) 
            + 0.0000528383579699300023*(taut3) 
            + 1.39398969845072005e-9*(pit7)*(taut_inv5*taut_inv) 
            + taut_inv7*(2.02584984300584983e-6*(pit3*pit)*taut_inv
            + 1.01874413933127995e-8*(pit7)*(taut_inv4) 
            + 0.000607063015658739955 - 0.000283190801238040004*taut_inv2)
            
            + pit14taut_inv17*(1.44400475720615078e-17*(pit3*pit3)*(taut_inv7*taut_inv5)
            - 3.33001080055983015e-19*pit7*pit*taut_inv7*taut_inv7
            + pit14taut_inv17*(3.58428679202129995e-22*(pit)*(taut_inv5)
            - 7.6373766822105502e-22*taut_inv4 
            - 5.65070932023524029e-23*pit2*taut_inv5*taut_inv 
            + 2.99318679335865594e-24*(pit3)*(taut_inv7))))

def iapws97_d2G_dpi2_region1(tau, pi):
    pit = 7.1 - pi
    taut = tau - 1.222
    taut_inv = 1.0/taut

    pit2 = pit*pit
    pit3 = pit2*pit
    pit6 = pit2*pit2*pit2
    pit15 = pit6*pit6*pit3
    
    taut2 = taut*taut
    taut5 = taut2*taut2*taut
    
    taut_inv2 = taut_inv*taut_inv
    taut_inv3 = taut_inv*taut_inv2
    taut_inv8 = taut_inv3*taut_inv3*taut_inv2
    taut_inv29 = taut_inv8*taut_inv8*taut_inv8*taut_inv3*taut_inv2
    return (-0.00001696247879118719984*pit - 0.0001900778690703240094*pit*taut_inv3*taut_inv - 5.112307687206180469e-9*pit*taut5*taut + 0.00009532278781397400275*taut - 0.0006000356158605199813 - 7.820546747472120262e-6*pit2*taut_inv2 - 0.00002691033828959999829*pit2*taut_inv3*taut_inv2 - 1.721007592550879906e-12*pit2*taut5*taut5 - 8.103399372023399331e-6*pit2*pit*taut_inv8 - 8.828369066169199338e-6*taut2*taut - 0.0009436864214653399542*taut_inv3 - 9.757927889155040732e-9*pit6*taut_inv3*taut_inv3 - 7.131208975318960007e-8*pit6*taut_inv8*taut_inv3 - 1.453899925951880029e-15*taut5*taut5*taut5*taut2 - 2.888009514412301439e-16*pit15*pit2*pit2*taut_inv29 + 7.326023761231625652e-18*(pit15*pit6)*(taut_inv29*taut_inv2) + 2.138465471018954057e-20*(pit15*pit6*pit6)*(taut_inv29*taut_inv8*taut_inv) - 1.039443169686176943e-20*pit15*pit6*pit6*pit*taut_inv29*taut_inv8*taut_inv2 + 1.695212796070572203e-21*(pit15*pit6*pit6*pit2)*(taut_inv29*taut_inv8*taut_inv3) - 9.278879059411833807e-23*pit15*pit15*taut_inv29*taut_inv8*taut_inv3*taut_inv)


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
    pi3 = pi*pi2
    pi4 = pi*pi3
    pi5 = pi2*pi3
    pi6 = pi3*pi3
    pi7 = pi6*pi
    pi9 = pi7*pi2
    pi15 = pi9*pi6
    pi19 = pi9*pi9*pi
    pi23 = pi9*pi9*pi5
    taut2 = taut*taut
    taut3 = taut*taut2
    taut4 = taut2*taut2
    taut6 = taut3*taut3
    taut7 = taut*taut6
    taut10 = taut4*taut6
    taut20 = taut10*taut10
    taut35 = taut20*taut10*taut3*taut2
    taut36 = taut*taut35
    taut50 = taut35*taut10*taut2*taut3
    return (pi19*(-2.93678005497663003e-14*pi3*taut35*taut4 + 0.0000832192847496054092*pi2*taut36*taut10*taut7 
            - 1.24017662339841913e-24*pi*taut10*taut10*taut) - 6.05920510335077989*pi15*pi2*taut50*taut7
            + pi4*taut7*(1.78287415218792009e-7*pi4*taut6 + 0.0000114610381688305001)
            - 0.000066065283340406*pi*taut - 0.000378979750326299998*pi*taut2
            - 0.00787855544867100029*pi*taut4 - 0.0875945913011459965*pi*taut7
            - 0.0000533490958281740028*pi*taut36 
            + taut10*(taut20*taut10*(taut10*(1.71088510070543998*pi15 
                            - 0.0000226487297378903988*taut4*taut4*pi23)
                      - 0.0000840049353964159951*taut4*taut4*pi19)
            # Cannot factor out taut here due to numerical issues
            + taut10*(1.32995316841867198e-15*taut20*pi23 - 1.29412653835175996e-9*taut7*taut2*pi15
                      + 1.78371690710842005e-23*pi19 + 1.75410265428146418e-27*taut6*pi23
                      - 0.27262789705017304*taut4*taut*pi6)
             
             - 0.012702883392812999*taut6*pi5
            - 1.00181793795109993e-8*taut4*pi9 - 8.83526622937069987e-6*taut*pi6
            - 1.02347470959289996e-12*pi9) + 9.00496908836719986e-11*taut4*taut4*pi7
            + 1.31612001853305008e-6*taut*pi2 - 3.15389238237468004e-9*taut*pi3
            - 0.0178348622923579989*taut - 0.0000968330317157100001*pi2*taut3
            - 0.00451017736264439952*pi2*taut6 - 0.122004760687946995*pi2*taut35
            + 6.14452130769269999e-8*pi2 - 4.13416950269890026e-17*pi6
            - 1.00288598706366e-10*pi5*taut3 - 143.374451604623999*pi5*taut35
            - 65.8490727183984035*pi7*taut36 + 1.04069652101739995e-18*pi9*taut4
            + 6.1258633752463995e-12*pi19*taut35 + 5.11628714091400033e-8*taut2*pi3
            - 0.0459960136963650026*taut2 + 1.92901490874028006e-6*taut3*pi3
            - 0.0575812590834320001*taut3 - 0.0503252787279300021*taut6 - 0.00177317424732129992)

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

### Region 3 subregion density calls


#print (_boundary_3rx(22E6), 0)


REGION_3A = 1
REGION_3B = 2
REGION_3C = 3
REGION_3D = 4
REGION_3E = 5
REGION_3F = 6
REGION_3G = 7
REGION_3H = 8
REGION_3I = 9
REGION_3J = 10
REGION_3K = 11
REGION_3L = 12
REGION_3M = 13
REGION_3N = 14
REGION_3O = 15
REGION_3P = 16
REGION_3Q = 17
REGION_3R = 18
REGION_3S = 19
REGION_3T = 20
REGION_3U = 21
REGION_3V = 22
REGION_3W = 23
REGION_3X = 24
REGION_3Y = 25
REGION_3Z = 26

REGION_3_AS_LETTERS = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E', 6: 'F', 7: 'G',
                       8: 'H', 9: 'I', 10: 'J', 11: 'K', 12: 'L', 13: 'M',
                       14: 'N', 15: 'O', 16: 'P', 17: 'Q', 18: 'R', 19: 'S',
                       20: 'T', 21: 'U', 22: 'V', 23: 'W', 24: 'X', 25: 'Y', 26: 'Z'}

def iapws97_region_3(T, P):
    r'''Identify the subregion in the IAPWS-97 region 3 given a temperature and 
    pressure point. No cheking that the original point is in region 3 is 
    performed.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]

    Returns
    -------
    subregion : int
        One of 1 to 26 inclusive, [-]

    Examples
    --------
    >>> iapws97_identify_region_TP(648.6, 22.5e6)
    3
    '''
    logP_MPa = log(P*1e-6)
    logP_MPa_inv = 1.0/logP_MPa
    region = 0
    if P > 100e6:
        return region
    Tab = iapws97_boundary_3ab(logP_MPa, logP_MPa_inv)
    if P > 40e6: #  and P <= 100e6
        if T <= Tab:
            return REGION_3A
        return REGION_3B
    Tcd = iapws97_boundary_3cd(P)
    Tef = iapws97_boundary_3ef(P)
    if P > 25e6 and P <= 40e6:
        if T <= Tcd:
            return REGION_3C
        elif Tcd < T <= Tab:
            return REGION_3D
        elif Tab < T <= Tef:
            return REGION_3E
        else:
            return REGION_3F
    
    Tgh = iapws97_boundary_3gh(P)
    Tij = iapws97_boundary_3ij(P)
    Tjk = iapws97_boundary_3jk(P)
    
    if P > 23.5E6 and P <= 25e6:
        if T <= Tcd:
            return REGION_3C
        elif Tcd < T and T <= Tgh:
            return REGION_3G
        elif Tgh < T and T <= Tef:
            return REGION_3H
        elif Tef < T and T <= Tij:
            return REGION_3I
        elif Tij < T and T<= Tjk:
            return REGION_3J
        return REGION_3K
    if 23e6 < P <= 23.5e6:
        if T <= Tcd:
            return REGION_3C
        elif Tcd < T <= Tgh:
            return REGION_3L
        elif Tgh < T <= Tef:
            return REGION_3H
        elif Tef < T <= Tij:
            return REGION_3I
        elif Tij < T <= Tjk:
            return REGION_3J
        return REGION_3K
    Tmn = iapws97_boundary_3mn(P)
    Top = iapws97_boundary_3op(logP_MPa, logP_MPa_inv)
    if 22.5e6 < P <= 23e6:
        if T <= Tcd:
            return REGION_3C
        elif Tcd < T <= Tgh:
            return REGION_3L
        elif Tgh < T <= Tmn:
            return REGION_3M
        elif Tmn < T <= Tef:
            return REGION_3N
        elif Tef < T <= Top:
            return REGION_3O
        elif Top < T <= Tij:
            return REGION_3P
        elif Tij < T <= Tjk:
            return REGION_3J
        else:
            return REGION_3K
    Trx = iapws97_boundary_3rx(P)
    Tqu = iapws97_boundary_3qu(P)
    Psat_643 = 21043367.318975247 # Psat_IAPWS(643.15)
    Tsat = Tsat_IAPWS(P)
    if Psat_643 < P <= 22.5e6:
        if T <= Tcd:
            return REGION_3C
        elif Tcd < T <= Tqu:
            return REGION_3Q
        elif Trx < T <= Tjk:
            return REGION_3R
        elif T > Tjk:
            return REGION_3K
        else:
            Tuv = iapws97_boundary_3uv(P)
            Twx = iapws97_boundary_3wx(logP_MPa, logP_MPa_inv)
            if 22.11e6 < P <= 22.5e6:
                if Tqu < T <= Tuv:
                    return REGION_3U
                elif Tuv < T <= Tef:
                    return REGION_3V
                elif Tef < T <= Twx:
                    return REGION_3W
                elif Twx < T <= Trx:
                    return REGION_3X
            if 22.064e6 < P <= 22.11e6:
                if Tqu < T <= Tuv:
                    return REGION_3U
                elif Tuv < T <= Tef:
                    return REGION_3Y
                elif Tef < T <= Twx:
                    return REGION_3Z
                elif Twx < T <= Trx:
                    return REGION_3X
            if T > Tsat:
                if Psat_643 < P <= 2.190096265E7:
                    #if T < Trx:
                    return REGION_3X
                if 2.190096265E7 < P <= 22.064e6:
                    if T <= Twx:
                        return REGION_3Z
                    elif Twx < T <= Trx:
                        return REGION_3X
            else:
                if Psat_643 < P <= 2.193161551E7:
                    return REGION_3U
                if 2.193161551E7 < P <= 22.064E6:
                    if Tqu < T <= Tuv:
                        return REGION_3U
                    return REGION_3Y
    if 20.5e6 < P <= Psat_643:
        if T <= Tcd:
            return REGION_3C
        elif Tcd < T <= Tsat:
            return REGION_3S
        elif Tsat <= T <= Tjk:
            return REGION_3R
        else:
            return REGION_3K
    P3cd = 19.00881189E6
    if P3cd < P <= 20.5e6:
        if T < Tcd:
            return REGION_3C
        elif Tcd < T < Tsat:
            return REGION_3S
        else:
            return REGION_3T
    Psat_623 = 16529164.25260448# Psat_IAPWS(623.15)
    if Psat_623 < P <= P3cd:
        if T <= Tsat:
            return REGION_3C
        return REGION_3T
    
        


def iapws97_region3_a(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*1e-08 - 0.085
    thetat = T*0.0013157894736842105 - 0.817
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.050226879086966297*pit + 0.040874841585674497*pit_inv + thetat*(0.47468639786331202*pit_inv + thetat*(-0.36964530819337699*pit + 1.1864681499791501*pit_inv + thetat*(1.3582570312914*(pit_inv2) + (thetat2)*(thetat*((thetat2)*((thetat2)*((thetat2)*(234105.65413187601*(pit_inv10) - 76705.194838085197*pit_inv10*pit_inv2) + 28277.661724328598*(pit_inv5) - 26989.395617661299*pit_inv5*pit_inv3 + 6280.0804934568896*(pit_inv10) + 572.61674081061597*(pit_inv10*pit_inv2)) - 2424.3152002952302*pit_inv2*pit_inv2 - 156.237904341963*pit_inv5*pit_inv3) + 44.272952105831401*(pit_inv3)) + 26.698704085604*(pit_inv5) + 0.216867826045856*(pit_inv5*pit_inv3) - 0.0253321069529674*pit_inv10 + 0.00110879558823853*(pit_inv10*pit_inv2)) + 1.7935760401998899*(pit_inv3)) + 0.079744179390101699*(pit2) + 0.45318626168577397*(pit_inv2)) + 0.19526677045264301 - 0.0122494831387441*pit_inv3 + 0.0011673222766826099*(pit_inv5) - 0.00018040710008550501*pit_inv5*pit_inv) + 0.54698726572754897 + 0.0063382803752842004*(pit2) - 0.0059322348901834198*pit_inv2 + 0.00043521732302273302*(pit_inv3))
    V = main*0.0024
    return 1.0/V
    

def iapws97_region3_b(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*1e-08 - 0.28
    thetat = T*0.0011627906976744186 - 0.779
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit3 = pit*pit2
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.102845919373532*pit + 0.0122261479925384*pit_inv + thetat*(thetat*(-0.49267663758928398*pit + 2.1635605769293802*pit_inv + thetat*(thetat*(thetat*(thetat*((thetat2)*((thetat2)*((thetat2)*(-29103.2084950276*pit_inv10*thetat2 + 41.688712601056501*(pit_inv10*pit_inv2)) - 1406.99677420738*pit_inv5 - 0.082767047000362096*pit_inv10*pit_inv2) + 361.18245261214901*(pit_inv5) + 140.24499760965799*(pit_inv5*pit_inv) - 111.422582236948*pit_inv5*pit_inv3 + 0.048365198219705897*(pit_inv10)) + 294.00250933851498*(pit_inv5*pit_inv)) + 23.960066025616101*(pit_inv2) - 50.667329572163702*pit_inv3 - 4.2559780405863199*pit_inv2*pit_inv2 - 344.38415881145897*pit_inv5 - 0.020230008390401399*pit_inv5*pit_inv) + 171.346792457471*(pit_inv2*pit_inv2)) - 41.375495701104199*pit_inv3) + 6.0881736840178498*(pit_inv2) - 0.24046253507852999*pit3 - 0.041637529016623598*pit_inv3 - 0.00202023902676481*pit_inv2*pit_inv2) - 0.116892827834085 + 0.00151140509678925*(pit_inv3) + 0.128369435967012*(pit2*pit2)) + 0.398198903368642 + 0.065554045640679001*(pit2) - 0.00057221296556902296*pit_inv2 + 6.9134608500033399e-6*(pit_inv3) - 0.026979818031007501*pit2*pit2)
    V = main*0.0041
    return 1.0/V
    

def iapws97_region3_c(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*2.5e-08 - 0.259
    thetat = T*0.0014492753623188406 - 0.903
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.156531975531713*pit + 0.0158646812591361*pit_inv + thetat*(0.70790633624184296*pit_inv + thetat*(11.770743004815801*pit + 12.601622514657*pit_inv + thetat*(thetat*(thetat*(thetat*(thetat*(thetat*((thetat2)*(-79389204.982125103*pit_inv10 + 32258310.340326902*(pit_inv10*pit_inv2)) - 3820310.2057081298*pit_inv2*pit_inv2 + 7912143.6522279195*(pit_inv5*pit_inv) - 899732.52990737697*pit_inv10 + 27671.3458847564*(pit_inv10*pit_inv2)) + 1232904.23502494*(pit2) - 1070777.1666086901*pit3 - 833426.56321285095*pit_inv5 + 175336.67532249901*(pit_inv5*pit_inv3)) + 2297.8474234507198*(pit_inv5*pit_inv3) - 342.41606509536302*pit_inv10 + 3.1196778876303002*(pit_inv10*pit_inv2)) + 3775.1566896695099*(pit_inv2) + 95.319300321738794*(pit_inv5*pit_inv3)) + 234.604891591616*(pit_inv2) - 65.950886355576699*pit_inv5) - 44.017020394964497*pit2 + 31.032749849200801*(pit_inv3)) - 17.810058818913699 + 0.064573468058329198*(pit_inv2*pit_inv2)) + 0.67654426899910103 - 0.18644246747194901*pit2 + 3.1993334584420903e-5*(pit_inv5) + 0.0438319858566475*(pit5*pit3)) + 0.73614365577215202 + 0.084014365386044704*(pit2) - 0.00089299671848372397*pit_inv2 - 0.0240650039730845*pit3 + 4.0639884847007902e-5*(pit_inv3))
    V = main*0.0022
    return 1.0/V
    

def iapws97_region3_d(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*2.5e-08 - 0.559
    thetat = T*0.0014492753623188406 - 0.939
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit3 = pit*pit2
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.032665012142638297*pit + 0.0050083391537212099*pit_inv + thetat*(0.38784248299841101*pit_inv + thetat*(thetat*(thetat*(thetat*(-1385.35367777182*pit_inv + thetat*(4980.4417172787698*pit + thetat*(thetat*((thetat2)*((thetat2)*((thetat2)*(-41125421794.6539*pit_inv10 + 1153711331204.97*(pit_inv10*pit_inv2)*(thetat2)) - 12712303.684593201*pit_inv10*pit_inv2) - 19541952.506071299*pit_inv5*pit_inv3 + 1443694.8990905299*(pit_inv10) + 508.05887480834502*(pit_inv10*pit_inv2)) + 2240407.5442698798*(pit_inv5*pit_inv) - 68931.508793315807*pit_inv5*pit_inv3 - 20.3578994462286*pit_inv10) + 125031.83535173599*(pit_inv2*pit_inv2) - 385294.21355528902*pit_inv5 - 22.177428114603799*pit_inv5*pit_inv3 - 0.0021499135204754499*pit_inv10*pit_inv2) + 3163.7351056401499*(pit_inv5*pit_inv) + 0.00277211346836625*(pit_inv10) + 3.1521038953880098e-5*(pit_inv10*pit_inv2)) - 348.15320341466298*pit_inv5) + 225.66051751243799*(pit_inv3) - 1.56481703640525e-6*pit_inv10 - 4.5248484717164498e-10*pit_inv10*pit_inv2) + 6.23449786243773e-6*(pit_inv5*pit_inv3)) + 1.7194625206874199 + 0.096812367845584099*(pit_inv3) - 0.00040421385283399602*pit_inv5 + 2.4155480603397201e-11*(pit_inv10)) - 0.029962841081922899*pit_inv2 + 0.00013464838327108901*(pit_inv2*pit_inv2) - 4.3670134792235601e-6*pit_inv5) + 0.87074524597177305 - 0.000190102435341872*pit_inv2 + 0.0055147802276508699*(pit3) + 1.3520370009940299e-7*(pit_inv2*pit_inv2) - 1.97805728776273e-16*pit_inv10)
    V = main*main
    V *= V*0.0029
    return 1.0/V
    

def iapws97_region3_e(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*2.5e-08 - 0.587
    thetat = T*0.0014084507042253522 - 0.918
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.54467412487890998*pit - 0.015749651617430801*pit_inv + thetat*(thetat*(thetat*(thetat*(2045.29931318843*pit - 1043.9079421301101*pit_inv2 + (thetat2)*(-22834.235932875199*pit + thetat*(thetat*((thetat2)*((thetat2*thetat2)*((thetat2)*(79497740233.560303*(pit_inv10) - 114328360753.44901*pit_inv10*pit_inv2) + 5353641749.6012697*(pit_inv10) + 715815808.40472102*(pit_inv10*pit_inv2)) - 1117963.81424162*pit_inv5*pit_inv3 + 665695.90883625206*(pit_inv10)) - 142586.073991215*pit_inv5*pit_inv3) - 266627.75039034098*pit_inv3 + 92.223056342143707*(pit_inv5*pit_inv3)) + 47599.266771712399*(pit_inv3) - 6699.8923907049102*pit_inv5 + 8961.2162964075997*(pit_inv5*pit_inv) - 9.0398366869115705e-5*pit_inv10) - 33.973132597771297*pit_inv2*pit_inv2) + 123.654999499486*(pit_inv2) + 3.7653100201571999e-12*(pit_inv10)) - 34.193183591040501*pit2 - 1.2052311155227799*pit_inv3 + 0.0045124253848683399*(pit_inv2*pit_inv2)) + 1.7837346287390301 + 0.30563840482826499*(pit_inv2)) + 0.685331118940253 + 0.41319748151589902*(pit2) - 0.00015331495438652401*pit_inv2)
    V = main*0.0032
    return 1.0/V
    

def iapws97_region3_f(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = sqrt(P*2.5e-08 - 0.587)
    thetat = T*0.0013698630136986301 - 0.891
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (-0.00141987303638727*pit*thetat_inv + thetat*(2.69251915156554*pit + thetat*(-30.020869577178299*pit*thetat + 34.974181585872202*pit - 16.517557195908601) + 2.1410775923648599 - 8.3909127728616895*pit2 - 1.4977453386065001*pit5 + 1.5120553127513301*(pit5*pit2)) - 0.00100615977450049*thetat_inv + 0.99996914025219197 - 1.31546288252539*pit2 + 6.0130719366876297e-6*(thetat_inv2) + 1.52115067087106*(pit3) - 0.00059109920647890904*pit3*thetat_inv2 + 1.81545608337015e-10*(pit3)*(thetat_inv5) - 2.5175654779232502e-8*thetat_inv3 + 2.5295647066322501e-5*(pit2*pit2)*(thetat_inv3) + 1.0072626520378601e-15*(pit5)*(thetat_inv5*thetat_inv3) - 7.9394097056296902e-10*pit5*pit*thetat_inv5*thetat_inv - 0.000150290891264717*pit5*pit2*thetat_inv2*thetat_inv2 + 4.7094260622165203e-6*(pit10)*(thetat_inv5*thetat_inv) + 0.00060437464020126502*(pit10*pit2)*(thetat_inv2*thetat_inv2) - 9.1162788626607693e-9*pit10*pit2*thetat_inv5*thetat_inv3 + 1.95049710391712e-13*(pit10*pit2)*(thetat_inv10) - 0.00091929673666610605*pit10*pit2*pit2*thetat_inv2*thetat_inv2 - 1.3779607079840901e-5*pit10*pit2*pit2*thetat_inv5*thetat_inv - 3.0306390804340398e-7*pit10*pit2*pit2*thetat_inv5*thetat_inv3 + 6.1091697358298098e-12*(pit10*pit2*pit2)*(thetat_inv10) - 2.2513293390013599e-16*pit10*pit2*pit2*thetat_inv10*thetat_inv2 + 7.5325947989869902e-7*(pit10*pit3*pit3)*(thetat_inv5*thetat_inv3) + 6.3928822313254498e-10*(pit10*pit3*pit3)*(thetat_inv10) + 7.5614029435161407e-9*(pit10*pit5*pit3)*(thetat_inv10) - 4.0032147868292898e-13*pit10*pit5*pit3*thetat_inv10*thetat_inv2 + 2.6958601059187399e-5*(pit10*pit5*pit3*pit2)*(thetat_inv5*thetat_inv) - 2.3761238114053899e-8*pit10*pit5*pit3*pit2*thetat_inv10 - 9.1208205403489092e-12*pit10*pit5*pit3*pit2*thetat_inv10*thetat_inv2 - 7.3282813515783902e-11*pit10*pit10*pit2*thetat_inv10*thetat_inv2 - 0.00040573553273032199*pit10*pit10*pit2*pit2*thetat_inv2*thetat_inv2 + 2.4199557830666001e-10*(pit10*pit10*pit2*pit2)*(thetat_inv10*thetat_inv2) + 1.8942414349801099e-10*(pit10*pit10*pit5*pit3)*(thetat_inv10*thetat_inv2) - 4.8663296507456297e-10*pit10*pit10*pit10*pit2*thetat_inv10*thetat_inv2)
    V = main*main
    V *= V*0.0064
    return 1.0/V
    

def iapws97_region3_g(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4e-08 - 0.872
    thetat = T*0.0015151515151515152 - 0.971
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.11048023927260101*pit + thetat*(5.3651603187505899*pit + thetat*(30.433158444409301*pit_inv + thetat*(-2914.41872156205*pit - 91.078254013468097*pit_inv2 + (thetat2)*(thetat*(thetat*(thetat*(5932507979.59445*pit_inv + (thetat2)*((thetat2)*((thetat2)*(-7.1294938340821105e+18*pit_inv2 + (thetat2*thetat2)*(-3.6417406211079798e+27*pit_inv + (thetat2)*((thetat2)*((thetat2)*(-1.04578785289542e+36*pit_inv2 + 6.16338176535305e+39*(pit3) + 7.2537907205934803e+29*(pit_inv10) - 1.05549884548496e+28*pit_inv10*pit_inv2) + 8.0205671552837802e+31*(pit_inv2*pit_inv2) - 1.2088917586118e+38*pit5 - 2.8021431005410101e+30*pit_inv5*pit_inv + 4.962507048713e+24*(pit_inv10*pit_inv2)) + 6.1375422916861901e+27*(pit_inv5) - 9.2217276959610093e+22*pit_inv10) - 1.9578886571897101e+17*pit_inv10*pit_inv2) - 758642165988.27795*pit_inv10 + 9481808850.3208008*(pit_inv10*pit_inv2)) + 8.1839602452461196e+22*(pit5*pit) + 228646846221.83099*(pit_inv5*pit_inv3) - 1149872.3828058699*pit_inv10*pit_inv2) - 37954580.233648703*pit_inv5*pit_inv3) - 4997410.93010619*pit_inv5*pit_inv + 10755.5033344858*(pit_inv5*pit_inv3)) - 29861781.9828065*pit_inv3 + 1049154.0676958601*(pit_inv5) - 61.7718249205859*pit_inv5*pit_inv3 + 4.1220902065299602e-5*(pit_inv10*pit_inv2)) - 8375139317986550.0*pit10) + 135033.22728156499*(pit_inv2)) + 940781944.83582902*(pit5*pit3)) - 72.464414375850794) - 0.33769360965747097) + 0.92179140353246103 - 36727.966954544798*pit10)
    V = main*main
    V *= V*0.0027
    return 1.0/V
    

def iapws97_region3_h(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4e-08 - 0.898
    thetat = T*0.0015151515151515152 - 0.983
    pit_inv = 1.0/pit
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.21134640224085799*pit + 0.0014719965861807599*pit_inv + thetat*(thetat*(24.997175295749098*pit + 20.261848702557799*pit_inv + thetat*(thetat*((thetat2)*(thetat*(thetat*((thetat2)*((thetat2)*((thetat2)*(17195156812433.699*(pit_inv10) - 18546115498514500.0*pit_inv10*thetat2) + 7741354215.8708296*(pit_inv10*pit_inv2)) - 605971823.58500504*pit_inv10) - 6561744.2199959401*pit_inv5*pit_inv + 17768333.734819099*(pit_inv5*pit_inv3) + 1936.9655876492*(pit_inv10) + 0.056137967888757703*(pit_inv10*pit_inv2)) - 2120.1062070121998*pit_inv5*pit_inv3) - 234396.09169331301*pit_inv5*pit_inv - 170.875935679023*pit_inv5*pit_inv3 - 0.00143987128208183*pit_inv10) + 1394.99167345464*(pit_inv2*pit_inv2) + 13.5249306374858*(pit_inv5) + 11.017744362957499*(pit_inv5*pit_inv) + 1.11482975877938e-9*(pit_inv10)) - 2.1294625702140002*pit_inv5) - 0.15201104438964799*pit_inv3 + 0.17718916414581301*(pit_inv2*pit_inv2) + 1.5636221297739601e-5*(pit_inv5)) - 0.0070367093203638799*pit_inv3 - 3.9546432784610502e-14*pit_inv5*pit_inv3) + 0.89934551894423997 + 9.8191692299111303e-5*(pit_inv2) + 3.8785116807801003e-17*(pit_inv5*pit_inv3))
    V = main*main
    V *= V*0.0032
    return 1.0/V
    

def iapws97_region3_i(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = sqrt(P*4e-08 - 0.91)
    thetat = T*0.0015151515151515152 - 0.984
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    thetat5 = thetat3*thetat2
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (-0.0023530288573684901*pit*thetat_inv - 0.269226321968839*pit - 5.6662075717003204e-7*pit*thetat_inv2 - 4.4635205567874902e-12*pit*thetat_inv2*thetat_inv2 + thetat*(thetat*((thetat3)*((thetat5)*(259862256980408.0 - 2.2328627042235599e+21*pit5*thetat2 + 4.2273953705724101e+19*(pit5*pit3)) - 99226310037675.0*pit10*pit2*pit2) + 42227580030.4086*(pit10*pit5*pit3)) - 1.4862085792233299) - 2.3177966967562398*thetat_inv*pit2*pit2 + 4.8133713145289097*thetat_inv*(pit5) + 1.0690568435913601 + 9.22024992944392*(pit2) - 17.394256556222199*pit3 + 3.5763350550377198e-12*(pit3)*(thetat_inv5) - 0.00026705035107576803*pit2*pit2*thetat_inv2 + 7.0068178555622898e-6*(pit2*pit2)*(thetat_inv3) - 7.5353304697975201e-13*pit5*thetat_inv5*thetat_inv + 0.0064641293413649596*(pit5*pit2)*(thetat_inv3) - 1.18746004987383e-5*pit5*pit2*thetat_inv2*thetat_inv2 - 4.1058853633093698e-10*pit5*pit3*thetat_inv5*thetat_inv + 3.1369818047381199e-13*(pit10)*(thetat_inv5*thetat_inv3) - 0.013526863990502101*pit10*pit2*thetat_inv2*thetat_inv2 - 3.39823323754373e-6*pit10*pit2*thetat_inv5*thetat_inv + 1.6439533434503999e-24*(pit10*pit2)*(thetat_inv10*thetat_inv2) - 0.046395953375238497*pit10*pit2*pit2*thetat_inv2*thetat_inv2 + 1.84386437538366e-9*(pit10*pit2*pit2)*(thetat_inv5*thetat_inv3) - 7.2325251421162496e-15*pit10*pit2*pit2*thetat_inv10 + 0.00345570606200257*(pit10*pit5*pit3)*(thetat_inv5*thetat_inv) - 5.4084301862408299e-8*pit10*pit5*pit3*thetat_inv5*thetat_inv3 - 2.2262099845219699e-11*pit10*pit5*pit3*thetat_inv10 + 6.8816915443933499e-17*(pit10*pit5*pit3)*(thetat_inv10*thetat_inv2) + 9.2723798515367901e-10*(pit10*pit5*pit3*pit2)*(thetat_inv10) - 1.2697447877048701e-15*pit10*pit5*pit3*pit2*thetat_inv10*thetat_inv2 + 6.1267081201648904e-14*(pit10*pit10*pit2)*(thetat_inv10*thetat_inv2) - 0.00038366950263682198*pit10*pit10*pit2*pit2*thetat_inv5*thetat_inv3 - 7.22693924063497e-12*pit10*pit10*pit2*pit2*thetat_inv10*thetat_inv2 - 93197.689751108599*pit10*pit10*pit10*pit2*thetat_inv5 + 0.000374684572410204*(pit10*pit10*pit10*pit2)*(thetat_inv10) + 65.811054675947403*(pit10*pit10*pit10*pit3*pit3)*(thetat_inv5*thetat_inv3) - 0.0247690616026922*pit10*pit10*pit10*pit3*pit3*thetat_inv10)
    V = main*main
    V *= V*0.0041
    return 1.0/V
    

def iapws97_region3_j(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = sqrt(P*4e-08 - 0.875)
    thetat = T*0.0014925373134328358 - 0.964
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (-0.000728541958464774*pit*thetat_inv + 1.79058760078792e-6*pit*(thetat_inv2) + thetat*(-18.757613337170401*pit + thetat*(6223.2297178647304*thetat*(pit5*pit) - 198.704578406823*pit2*pit2) + 5.3061558192897902 + 24.357475537729002*(pit2)) - 0.00011137131739554*thetat_inv + 0.0019906087407184901*thetat_inv*(pit2) + 1.0034289242368499 - 0.000177040785499444*pit3*thetat_inv2 - 0.0025968038522712999*pit2*pit2*thetat_inv2 - 1.6102312131433301*pit5 - 0.0023626469284413801*pit5*thetat_inv2 + 7.3862779022428697e-5*(pit5)*(thetat_inv3) - 9.6075411670166893e-9*pit10*thetat_inv5*thetat_inv + 0.0076737378140421097*(pit10*pit2)*(thetat_inv3) - 5.1057226972048803e-11*pit10*pit2*thetat_inv5*thetat_inv3 + 1.46564542926508e-5*(pit10*pit2*pit2)*(thetat_inv5) - 7.1759073552674496e-10*pit10*pit2*pit2*thetat_inv5*thetat_inv3 + 6.6385546948525403e-15*(pit10*pit2*pit2)*(thetat_inv10) + 3.0902947427701301e-12*(pit10*pit3*pit3)*(thetat_inv10) - 4.6421630097170797e-16*pit10*pit5*pit3*thetat_inv10*thetat_inv2 - 2.36716126781431e-10*pit10*pit5*pit3*pit2*thetat_inv10 - 3.9049963796116099e-14*pit10*pit5*pit3*pit2*thetat_inv10*thetat_inv2 - 0.0042227178748249702*pit10*pit10*pit2*pit2*thetat_inv5*thetat_inv + 4.54652854268717e-12*(pit10*pit10*pit2*pit2)*(thetat_inv10*thetat_inv2) + 2.7092900272022802*(pit10*pit10*pit5*pit3)*(thetat_inv5) + 2.8391174235470601e-11*(pit10*pit10*pit5*pit3)*(thetat_inv10*thetat_inv2))
    V = main*main
    V *= V*0.0054
    return 1.0/V
    

def iapws97_region3_k(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4e-08 - 0.802
    thetat = T*0.0014705882352941176 - 0.935
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    thetat2 = thetat*thetat
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (-3.4470960548668601*pit - 0.00027761760697574801*pit*thetat_inv2 + 1.7034307284185001e-6*pit*(thetat_inv3) + 3.9472147136367802e-15*pit_inv*(thetat_inv5) + thetat*(22.133386244709499*pit + thetat*(-194.64611003707901*pit + thetat*(thetat*(3289.13873658481*(pit2) + (thetat2)*(37262.996737414702*pit_inv + (thetat2*thetat2)*(-401215699.57609898*pit_inv2 + (thetat2)*(48450147831.840599*(pit_inv2) - 29102685116444.398*thetat2)))) + 589.70277127742895) - 104.52963483027899) + 12.243316265660001) - 0.000879148916140706*thetat_inv + 0.84431786384433105 + 2.5583029857902702*(pit2) - 0.00181057560300994*pit2*thetat_inv2 - 6.96664158132412e-6*pit2*thetat_inv3 - 1.8084520914547001e-11*pit2*thetat_inv5*thetat_inv + 8.0835463977282498e-16*(pit2)*(thetat_inv5*thetat_inv3) + 4.7536162997023301e-7*(thetat_inv2) - 0.0039568892342125*pit5*thetat_inv3 - 6.6187679255803405e-7*pit5*thetat_inv5*thetat_inv - 1.73270241249904e-19*pit5*thetat_inv10*thetat_inv2 + 3.8371940902555598e-5*(pit5*pit)*(thetat_inv5) + 1.6075110746495799e-9*(pit5*pit)*(thetat_inv5*thetat_inv3) - 4.00879935920517e-14*pit5*pit*thetat_inv10 + 6.0420329981913199e-18*(pit5*pit)*(thetat_inv10*thetat_inv2) - 3.8043640701245203e-15*thetat_inv5*thetat_inv - 6.4956544670245702e-15*pit5*pit3*thetat_inv10*thetat_inv2 - 1.4909532850600001e-12*pit10*thetat_inv10*thetat_inv2 + 5.4144937732958098e-9*(pit10*pit2)*(thetat_inv10) - 3.6979437416866603e-30*thetat_inv10*thetat_inv2)
    V = main*0.0077
    return 1.0/V
    

def iapws97_region3_l(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.166666666666667e-08 - 0.908
    thetat = T*0.0015384615384615385 - 0.989
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    thetat5 = thetat3*thetat2
    thetat10 = thetat5*thetat5
    main = (-0.198068404154428*pit + 2.6583961888553e-5*pit_inv + thetat*(0.025339239288975399*pit_inv + thetat*(thetat*(-214.443041836579*pit_inv + thetat*(thetat*(thetat*(thetat*(thetat*((thetat2)*(-306367307532219.0*pit_inv2 + (thetat2)*(4.9423723717971798e+20 + (thetat2)*((thetat2)*(-1.4141534988114001e+30*pit + (thetat2)*((thetat2)*((thetat2)*((thetat2)*(5.21635864527315e+34*(pit_inv5*pit_inv3) - 4.87095672740742e+54*pit_inv5*pit_inv3*thetat10*thetat2 - 3.8145826048995503e+32*pit_inv10) + 4.1386518684890801e+26*(pit_inv10*pit_inv2)) - 6.9595362234882901e+32*pit_inv3 - 7.5896694638775806e+22*pit_inv10*pit_inv2) + 2.7070611108523801e+29*(pit_inv3) - 4.4435947874629502e+22*pit_inv5*pit_inv3 + 5.5492387028966697e+18*(pit_inv10*pit_inv2)) - 1.0810548079647099e+24*pit_inv2*pit_inv2 - 188277213604704.0*pit_inv10*pit_inv2) + 1.1666212121932199e+32*(pit5*pit) + 5294829964228630.0*(pit_inv5*pit_inv3) - 815038000738.06006*pit_inv10 + 2607020586.4753699*(pit_inv10*pit_inv2)) - 495017809506.71997*pit_inv5*pit_inv3 - 4.4570336919694501e+32*pit10) + 22609563.143717401*(pit_inv5*pit_inv3) + 6.4279493237369401e+32*(pit10*pit2*pit2)) - 714430.20993754698*pit_inv5*pit_inv) + 7774514.3796098996*(pit_inv2*pit_inv2)) - 0.0123239564600519*pit_inv5*pit_inv3) - 10.0752127917598*pit_inv5) + 0.127868634615495*(pit_inv5) - 3158749762715330.0*pit10) + 72.155916336135405*(pit_inv2) - 2.1285716942348398*pit_inv3) + 33.840122250919102 + 0.11060902747228001*(pit_inv2)) + 2.2318404310169999 - 99.386242161365104*pit2 - 3.5757858116965902e-6*pit_inv3 + 47313.790987276501*(pit5)) + 0.93784660148966703 + 125.07053414273101*(pit2*pit2) - 996.473529004439*pit5)
    V = main*main
    V *= V*0.0026
    return 1.0/V
    

def iapws97_region3_m(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 1.0
    thetat = sqrt(sqrt(T*0.0015384615384615385 - 0.997))
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    main = (0.81138436348184695 + (thetat2)*((thetat3)*(-81456.820934687203*pit + thetat*(458384.82859394897*pit + thetat*(thetat*((thetat2)*(-5475783138.9909697*pit2 + (thetat2)*((thetat2)*(-170451090076.38501*pit + 185135446.82833701 + (thetat2*thetat2)*(157890366037614.0*pit + (thetat2)*(-2025305097487740.0*pit + (thetat2)*(1.70215539458936e+17*(pit2) + (thetat2)*(-821698160721956.0 + (thetat2*thetat2)*(2.3341586947851002e+17 - 6.0007993458680299e+22*pit3 + 5.9458438227338403e+24*(pit2*pit2) + (thetat2*thetat2)*(1.88813911076809e+21*pit + (thetat2*thetat2)*(-3.2942192395146001e+21 - 1.37570282536696e+25*pit2 + 1.8150899630390199e+27*(pit3) - 3.4686512276835299e+29*pit2*pit2 - 2.1196114877426e+37*pit5*pit3 - 1.2861789988767499e+48*pit10*pit2*pit2 + 4.79817895699239e+64*(pit10*pit10*pit2*pit2)) + 1.11052244098768e+35*(pit5*pit3) + 2.9113395860250301e+45*(pit10*pit2*pit2)) + 1.8946127934949199e+39*(pit10*pit2) - 8.1009342884264494e+45*pit10*pit3*pit3) - 7.9526024187230606e+23*pit5) + 6.3923490991874096e+41*(pit10*pit3*pit3)) + 3.6819392618356999e+59*(pit10*pit10*pit5*pit3)))) + 1850072455632.3899*(pit3)) + 200725701112386.0*(pit5)) + 939454935735.56299*(pit2*pit2) + 2.6657285643293799e+27*(pit10*pit2*pit2)) + 45373580.000427298*(pit2)) - 38575400038384.797*pit5*pit) - 65977456.760287397*pit3 - 15286114865.930201*pit2*pit2 - 560165667510.44604*pit5) + 7.9553765761342698e+31*(pit10*pit5*pit3*pit2)) - 5681.99310990094*pit3 - 17865719817.2556*pit5*pit3)
    V = main*0.0028
    return 1.0/V
    

def iapws97_region3_n(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.976
    thetat = T*0.0015384615384615385 - 0.997
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat2 = thetat*thetat
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (-302.807107747776*pit + thetat*(232534.27270987601*pit + thetat*((-86987136466.276901*pit + thetat*(400849240129329.0*pit*thetat + 354542769185.67102))*(thetat2) - 792681.20713260002) + 1591.5874831459901) - 0.00089076330670130501*thetat_inv - 4402.0959940771399*thetat_inv*pit3 - 4.9311136203016203e-11*pit2*thetat_inv5 + 5.4127691156417598e-14*(pit2)*(thetat_inv5*thetat_inv) + 2.4056080832171301e-7*(thetat_inv2) - 0.0064306413263692502*pit3*thetat_inv3 + 7.0541210077369902e-12*(pit3)*(thetat_inv5*thetat_inv) - 3.34952758812999e-19*pit3*thetat_inv5*thetat_inv3 - 6.0724664397089301e-24*pit3*thetat_inv10 + 6.14869006573609e-31*(pit3)*(thetat_inv10*thetat_inv2) + 0.0022001990172961501*(pit2*pit2)*(thetat_inv2*thetat_inv2) - 1.5864969989454301e-6*pit2*pit2*thetat_inv5 + 2.5858588789748602e-9*(pit2*pit2)*(thetat_inv5*thetat_inv) + 5.8223866704894202e-28*(pit2*pit2)*(thetat_inv10*thetat_inv2) + 62.915414901504803*(pit5)*(thetat_inv3) - 4.0235211523449402e-19*pit5*thetat_inv10 + 135.14731861706099*(pit5*pit)*(thetat_inv3) - 7.44938506925544e-17*pit5*pit*thetat_inv10 + 3.9062836923846201e-23*(pit5*pit)*(thetat_inv10*thetat_inv2) - 0.52503742788609997*pit5*pit2*thetat_inv5 - 4.2153772609838902e-9*pit5*pit2*thetat_inv5*thetat_inv3 + 8.2144575825511904e-21*(pit5*pit2)*(thetat_inv10*thetat_inv2) + 1.8991720652623699e-13*(pit5*pit3)*(thetat_inv10) + 1.7727487236194601e-26*(thetat_inv5*thetat_inv3) + 4.0213796184277599e-15*(pit10)*(thetat_inv10*thetat_inv2) - 1.35031446451331e-32*thetat_inv10 - 0.039104816792964903*pit10*pit2*thetat_inv5*thetat_inv3 + 3.64975183508473e-6*(pit10*pit2)*(thetat_inv10) + 6.5171817187830098e-13*(pit10*pit2)*(thetat_inv10*thetat_inv2) + 2.80967799943151e-39*(thetat_inv10*thetat_inv2) - 2.1177335580305798e-8*pit10*pit2*pit2*thetat_inv10*thetat_inv2 + 0.0026495335438007201*(pit10*pit5*pit3)*(thetat_inv10*thetat_inv2))
    V = exp(main)*0.0031
    return 1.0/V
    

def iapws97_region3_o(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = sqrt(P*4.347826086956522e-08 - 0.974)
    thetat = T*0.0015384615384615385 - 0.996
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (0.0028907869214915001*thetat_inv + 0.244482731907223*thetat_inv*(pit2) + 1.38647388209306*thetat_inv*(pit2*pit2) + 1.4173349203098499e-24*(pit3)*(thetat_inv10) + 2.0137732541180301e-6*(pit2*pit2)*(thetat_inv2*thetat_inv2) - 5.8518840178277898e-9*pit2*pit2*thetat_inv5 - 5.9453920290143103e-18*pit2*pit2*thetat_inv5*thetat_inv3 - 3.5453385305947602e-29*pit2*pit2*thetat_inv10*thetat_inv2 - 7.35234770382342e-12*thetat_inv2*thetat_inv2 + 0.00137680878349369*(pit5)*(thetat_inv3) - 1.7395936508477201e-5*pit5*thetat_inv2*thetat_inv2 + 8.1489760580551298e-15*(pit5*pit)*(thetat_inv5*thetat_inv3) + 4.2559663135183898e-26*(pit5*pit2)*(thetat_inv10*thetat_inv2) - 0.0017184963895152099*pit5*pit3*thetat_inv2*thetat_inv2 + 1.3981474793024e-13*(pit5*pit3)*(thetat_inv5*thetat_inv3) - 3.8744911378775499e-18*pit5*pit3*thetat_inv10 + 1.1896057807201801e-11*(pit10)*(thetat_inv5*thetat_inv3) + 6.4189052951329603e-22*(pit10)*(thetat_inv10*thetat_inv2) + 1.2874602397971801e-35*(thetat_inv10*thetat_inv2) + 2.33907907347507e-8*(pit10*pit2*pit2)*(thetat_inv5*thetat_inv3) - 1.55282762571611e-18*pit10*pit2*pit2*thetat_inv10*thetat_inv2 + 3.7768264908914897e-9*(pit10*pit5*pit3*pit2)*(thetat_inv10) - 1.7409324776621299e-13*pit10*pit5*pit3*pit2*thetat_inv10*thetat_inv2 - 5.1672023657530198e-11*pit10*pit10*pit2*pit2*thetat_inv10*thetat_inv2)
    V = main*0.0034
    return 1.0/V
    

def iapws97_region3_p(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = sqrt(P*4.347826086956522e-08 - 0.972)
    thetat = T*0.0015384615384615385 - 0.997
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (thetat*(-1235.9234861013699*pit + 3246.64750281543*thetat + 116.033094095084) - 9.8282534201036601e-5*thetat_inv - 0.056140345001349498*thetat_inv*pit2 + 1.0514570085061199 + 236.31342539392401*(pit3) + 8.5667740164086898e-8*(pit3)*(thetat_inv3) + 0.0097250329235010896*(pit2*pit2)*(thetat_inv2) - 1.03001994531927*pit5*pit*thetat_inv2 - 2.15743778861592e-5*pit5*pit2*thetat_inv2*thetat_inv2 - 1.4965370619916199e-9*pit5*pit2*thetat_inv5 - 8.3445219829144506*pit5*pit3*thetat_inv2 + 0.586602660564988*(pit10)*(thetat_inv3) + 0.00294985697916798*(pit10*pit2)*(thetat_inv5) + 8.1625609594702108e-6*(pit10*pit2)*(thetat_inv5*thetat_inv) + 3.4348002210496797e-26*(pit10*pit2)*(thetat_inv10*thetat_inv2) + 10.776602703285301*(pit10*pit2*pit2)*(thetat_inv3) + 4.0095476380694099e-10*(pit10*pit2*pit2)*(thetat_inv5*thetat_inv3) + 7.1173046627658394e-17*(pit10*pit2*pit2)*(thetat_inv10) - 4.0944959913818202e-7*pit10*pit3*pit3*thetat_inv5*thetat_inv3 - 7.2912130775890196e-6*pit10*pit5*pit3*thetat_inv5*thetat_inv3 + 6.7710797093890902e-9*(pit10*pit5*pit3*pit2)*(thetat_inv10) + 6.0274597302297505e-8*(pit10*pit10*pit2)*(thetat_inv10) + 0.00179946628317437*(pit10*pit10*pit2*pit2)*(thetat_inv5*thetat_inv3) - 3.82323011855257e-11*pit10*pit10*pit2*pit2*thetat_inv10*thetat_inv2 - 0.000345042834640005*pit10*pit10*pit10*pit3*pit3*thetat_inv10*thetat_inv2)
    V = main*0.0041
    return 1.0/V
    

def iapws97_region3_q(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.848
    thetat = T*0.0015384615384615385 - 0.983
    pit_inv = 1.0/pit
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    main = (-0.0838165632204598*pit + 0.0010035565172151*pit_inv + thetat*(2.4779590841149202*pit + 0.33349145514351602*pit_inv + thetat*(1.0969757688887301*pit_inv + thetat*(-3191.1496900653301*pit + thetat*(thetat*(thetat*(thetat*(thetat*((thetat2)*(-1729857814.3333499*pit_inv10 - 82043.384325994994*pit_inv10*pit_inv2 + 47327151846.1586*(pit_inv10*pit_inv2)*(thetat2)) + 35176923.2729192*(pit_inv5*pit_inv3) - 3566.1702998249002*pit_inv10) + 32.860002543598*(pit_inv10)) - 775489.25998514402*pit_inv5*pit_inv - 0.080595002100541296*pit_inv10) + 99349.988382027397*(pit_inv5)) + 2256.8993916191798*(pit_inv2) - 6128.4281682008304*pit_inv2*pit_inv2) + 232.808472983776*(pit_inv3) - 0.64209417190456997*pit_inv2*pit_inv2) - 4.2857722747561402*pit_inv2 + 7.1034669196601803e-5*(pit_inv5)) - 0.0064359606067845602*pit_inv2) + 0.96191737937645205 - 1.42808220416837e-5*pit_inv2)
    V = main*main
    V *= V*0.0022
    return 1.0/V
    

def iapws97_region3_r(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.874
    thetat = T*0.0015384615384615385 - 0.982
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    thetat5 = thetat3*thetat2
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (thetat*(3.0530889006508901 + (thetat2)*(thetat*(thetat*(thetat*((thetat2)*(490112654.15421098*(pit_inv3) - 7014385996282.5801*pit_inv5*pit_inv3*thetat5*thetat) + 0.0014416595566086299*(pit_inv5*pit_inv3)) - 3997452.76971264 - 10433.403065402101*pit_inv3) + 393.09721470624498*(pit_inv3)) + 0.26197513536810901*(pit_inv3))) - 0.000147104222772069*thetat_inv + 1.0360274804340801 - 0.046492350440777798*pit3*thetat_inv2 + 5.6923371959374999e-12*(pit3)*(thetat_inv5*thetat_inv) - 8.30946716459219e-17*pit_inv3*thetat_inv3 + 0.015953672241120199*(pit5*pit3)*(thetat_inv5) - 5.3647956020181104e-7*pit5*pit3*thetat_inv5*thetat_inv3 + 3.9998879569316202e-13*(pit5*pit3)*(thetat_inv10) - 5.3540039651290598e-18*pit5*pit3*thetat_inv10*thetat_inv2 + 150764.97412551101*(pit10)*(thetat_inv2) - 14336.5406393758*pit10*thetat_inv3 + 546.49132352849097*(pit10)*(thetat_inv2*thetat_inv2) - 9.9345695784500592*pit10*thetat_inv5 + 0.066351314422445407*(pit10)*(thetat_inv5*thetat_inv) - 9.8343063671645403e-6*pit10*thetat_inv5*thetat_inv3 + 2.4424745385850599e-8*(pit10)*(thetat_inv10) + 2.70303248860217e-15*(pit10)*(thetat_inv10*thetat_inv2) - 3.37209709340105e-10*pit10*pit2*thetat_inv10*thetat_inv2 + 3.7750198002546901e-9*(pit10*pit2*pit2)*(thetat_inv10*thetat_inv2))
    V = main*0.0054
    return 1.0/V
    

def iapws97_region3_s(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.761904761904762e-08 - 0.886
    thetat = T*0.0015625 - 0.99
    pit_inv = 1.0/pit
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    thetat5 = thetat3*thetat2
    main = (-0.198339358557937*pit + thetat*(0.0012139997999321701*pit_inv + thetat*(1.88317043049455*pit_inv + thetat*(-1670.7350396206*pit_inv + thetat*(-65391.5627346115 + 18751.4491833092*(pit_inv2) + (thetat2)*((thetat2)*(-5041607241.3259001*pit_inv3 + 88458547.259613395*(pit_inv5) + (thetat5*thetat)*((thetat2)*((thetat2)*((thetat2)*((thetat2)*((thetat2)*(-6.1655261113579205e+45*pit2*pit2 + (thetat2*thetat2)*(6.0401220016344396e+49 + (thetat2*thetat2)*(-1.7598409016350101e+57*pit - 1.8566232754532401e+53*pit_inv2*pit_inv2 + 2.02281884477061e+58*(pit_inv5*pit_inv)*(thetat2*thetat2))) + 1.00415480000824e+31*(pit_inv10*pit_inv2) + 9.5089817042504195e+53*(pit10*pit2*pit2)) - 1.9154000182136699e+29*pit_inv10) - 5.3246661214025397e+22*pit_inv10*pit_inv2) + 4.37796099975134e+33*(pit2*pit2)) + 1.66540181638363e+22*(pit_inv5)) + 10561837780884700.0*(pit_inv5*pit_inv3))) - 313563.19766911102*pit_inv2*pit_inv2) + 1935687689.1779699*(pit5)) - 0.0624942093918942*pit_inv3 - 10917404.498782899*pit2*pit2) + 45621.341533807099*(pit3)) + 2.94885696802488 - 575.99125514438401*pit3) + 0.96596165059977501 + 3.5631488140398702*(pit3))
    V = main*main
    V *= V*0.0022
    return 1.0/V
    

def iapws97_region3_t(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*5e-08 - 0.803
    thetat = T*0.0015384615384615385 - 1.02
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat2 = thetat*thetat
    main = (-2.9100291578376098*pit + thetat*(6.6423511500903096 + (thetat2)*(thetat*(-2893.6623672720998 + (thetat2)*(thetat*(thetat*((thetat2)*(-829088246858.08301*pit + (thetat2)*(-3859232023098.48 + (thetat2)*(1.6046460868783398e+17*(pit2) + (thetat2)*((thetat2)*((thetat2)*((thetat2)*((thetat2)*((thetat2*thetat2)*((thetat2*thetat2)*((thetat2*thetat2)*(-3.4155204086064402e+50*pit5*pit2 - 7.0577262332637399e+64*pit10*pit10*pit2 - 4.4422736775830395e+71*pit10*pit10*pit10*pit2 - 2.81396013562745e+76*pit10*pit10*pit10*pit3*pit3) + 5.6902145441327e+57*(pit10*pit5*pit3*pit2) + 4.2843233862067799e+68*(pit10*pit10*pit10*pit2)) - 3.00475129680486e+60*pit10*pit10*pit5*pit3) + 1.6686117620014801e+52*(pit10*pit10*pit2*pit2)) - 6.5647528033941103e+35*pit10 - 7.0058454643311301e+47*pit10*pit10*pit2 - 6.6848129519680801e+50*pit10*pit10*pit10*pit2) - 3.2791059208652301e+30*pit5*pit2) + 3.5528604551230099e+38*(pit10*pit5*pit3)) + 3.58958955867578e+28*(pit10)) - 1.68776617209269e+26*pit10) + 2.4537564093705501e+23*(pit10)) - 2297462376236920.0*pit2*pit2 - 5.2725133970904698e+20*pit10) + 1566374275417.29*(pit3)) - 67707383068734.898*pit5*pit2) - 534686695.71346903*pit2) + 1105544467.9054301*(pit5*pit2)) + 196435.366560186*(pit3) + 38565900.164800599*(pit5*pit2))) + 1.5528724958626801 + 1.76814899675218*(pit2) - 1.78154560260006*pit2*pit2)
    V = main*0.0088
    return 1.0/V
    

def iapws97_region3_u(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.902
    thetat = T*0.0015384615384615385 - 0.988
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (0.00105581745346187*pit*(thetat_inv2) - 0.00012270822923564101*pit_inv*thetat_inv + thetat*(-106.201671767107*pit_inv + thetat*(thetat*(thetat*(thetat*(thetat*((thetat2)*((thetat2)*(-1.6011681327467599e+24*pit2 + (thetat2)*(9.03443213959313e+24*pit_inv + (thetat2)*(-6.9399627037085199e+27*pit_inv - 3.1443257755155201e+21*pit_inv5*pit_inv3 + 2.5992951084949901e+19*(pit_inv10) + 1.2208834925835501e+17*(pit_inv10*pit_inv2)) - 8.4340592684641799e+20*pit_inv5 + 1.5907964819684901e+20*(pit_inv5*pit_inv) - 8.7847358505008499e+17*pit_inv5*pit_inv3 - 8826669315646520.0*pit_inv10) + 222612779142211.0*(pit_inv5*pit_inv3) + 1042164686.08488*(pit_inv10)) + 8843876513378.3594*(pit_inv5) - 2169349169962.8501*pit_inv5*pit_inv + 1.04674840020929e+26*(pit5*pit3)) - 7.8175450769884603e+27*pit10*pit2*pit2) - 651903203602581.0*pit2) - 339.56761730342299*pit_inv5 + 2.2614596374788099e+21*(pit10*pit2)) + 276378438378930.0*(pit5)) + 11.4178193518022*(pit_inv3) + 677143292290.14404*(pit5) - 30142694798017.102*pit5*pit) + 7189.5756712785096) - 5.1025429423783704e-9*pit3*thetat_inv5 + 6.4891671896557497e-9*(thetat_inv3) - 0.152355388953402*pit5*thetat_inv2*thetat_inv2 + 0.0116862983141686*(pit5*pit)*(thetat_inv5) + 1.6971981388484e-8*(pit5*pit3)*(thetat_inv5*thetat_inv3) - 10801.690456013999*pit10*thetat_inv2*thetat_inv2 + 5361164.8360273801*(pit10*pit2)*(thetat_inv2*thetat_inv2) - 9.9062360193429495e-13*pit10*pit2*thetat_inv10*thetat_inv2 - 22770.046464391999*pit10*pit2*pit2*thetat_inv5*thetat_inv + 1.5100154888067e-5*(pit10*pit2*pit2)*(thetat_inv10) - 4.8873156577621002e-10*pit10*pit2*pit2*thetat_inv10*thetat_inv2)
    V = main*0.0026
    return 1.0/V
    

def iapws97_region3_v(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.96
    thetat = T*0.0015384615384615385 - 0.995
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (2.7603260114515101e-29*pit*(thetat_inv10) - 1.11526741826431e-35*pit*thetat_inv10*thetat_inv2 + 0.018458726111483699*pit_inv - 1.8821488234144798e-9*pit_inv*thetat_inv2 + thetat*(thetat*(thetat*(thetat*(thetat*(thetat*(-72368188562634800.0 + (thetat2)*((thetat2)*(-2.2344919405412401e+26 + (thetat2)*(7.4270572330273797e+26*(pit_inv3) - 1.03977184454767e+28*pit_inv5*thetat2) - 1.92359972440634e+22*pit_inv3 - 4.6813835890873197e+31*pit2*pit2 + 5.8779310562074801e+20*(pit_inv2*pit_inv2) - 6.9759575034739098e+18*pit_inv5 + 31308029991594400.0*(pit_inv5*pit_inv)) + 5131174628650.4404*(pit_inv5) - 62418400710.315804*pit_inv5*pit_inv) + 654144.37374993705*(pit_inv5) + 59461.976619346002*(pit_inv5*pit_inv)) - 25.912373638026899*pit_inv5*pit_inv) + 8206120.4864546899*(pit_inv2)) + 134856491567853.0*(pit3) + 51065511977436000.0*(pit2*pit2)) - 51.742968245060503*pit_inv2 - 7606674911832790.0*pit5 - 1.92824336984852e-6*pit_inv5) + 1.05006446192036e-9*(pit_inv5) + 2.4776139232905799e+26*(pit10*pit2*pit2)) - 1.3583040778266301e-6*thetat_inv2 + 2.80375725094731e-18*(pit_inv3)*(thetat_inv3) + 6.5244029334585999e-10*(pit2*pit2)*(thetat_inv5*thetat_inv) + 9.2699003653063902e-30*(pit_inv2*pit_inv2)*(thetat_inv5*thetat_inv) - 4.3667703405165502e-42*pit_inv2*pit_inv2*thetat_inv10 + 1.1956313554066601e-48*(pit_inv2*pit_inv2)*(thetat_inv10*thetat_inv2) + 3.5925221360411398e-26*(pit_inv5*pit_inv)*(thetat_inv3) - 3.5707866820337699e-55*pit_inv5*pit_inv*thetat_inv10*thetat_inv2 - 4.1724798698682099e-19*pit5*pit3*thetat_inv10*thetat_inv2 + 1.7744174292404301e-61*(pit_inv5*pit_inv3)*(thetat_inv10*thetat_inv2) + 31254567775610.398*(pit10)*(thetat_inv2) - 4.15652812061591e-55*pit_inv10*thetat_inv5*thetat_inv3 - 100375333864186.0*pit10*pit2*thetat_inv3)
    V = main*0.0031
    return 1.0/V
    

def iapws97_region3_w(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.959
    thetat = T*0.0015384615384615385 - 0.995
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    pit_inv10 = pit_inv5*pit_inv5
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    thetat5 = thetat3*thetat2
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (0.92326135790147001*pit*thetat_inv + 2.7170023573989301e-15*pit_inv*(thetat_inv2*thetat_inv2) + 2.3741673261664401e-27*pit_inv*(thetat_inv5*thetat_inv3) + thetat*(-90.788621348359996*pit_inv + thetat*(thetat*(-1033.08436323771*pit_inv3 + (thetat3)*((thetat2)*(-96110924.098574698*pit_inv5*pit_inv + (thetat5*thetat)*(-1.5854860965500201e+18*pit_inv5*pit_inv3 - 89446035500.552597*pit_inv10*pit_inv2) + 22827.6853990249*(pit_inv5*pit_inv3) + 0.109892402329239*(pit_inv10) - 5.8621913381701603e-8*pit_inv10*pit_inv2) - 0.0575368389425212*pit_inv5*pit_inv3) + 0.725937724828145*(pit_inv2*pit_inv2)) + 3219887.6763638901*(pit2) + 579.51404176570998*(pit_inv2) + 6.15762068640611e-9*(pit_inv5*pit_inv)) + 156.792067854621 - 0.066255281634216803*pit_inv2) - 5.9786598842257703*thetat_inv*pit2 - 4.7110372549807704e-13*thetat_inv*pit_inv2*pit_inv2 + 5.3116803751977397e-31*thetat_inv*(pit_inv10) + 4.9342908604698102e-8*(pit3)*(thetat_inv5) - 3.9944139004220299e-30*pit3*thetat_inv10*thetat_inv2 + 1.8776852576368201e-39*(pit_inv3)*(thetat_inv10) - 3.4082129141971899e-7*pit5*thetat_inv5*thetat_inv - 2.0761028465413701e-12*pit5*thetat_inv5*thetat_inv3 + 8.1203698337056496e-20*(pit5)*(thetat_inv10) - 4.0627428665262499e-45*pit_inv5*thetat_inv10 - 6.3498798119066896e-25*pit_inv5*pit_inv*thetat_inv3 + 3.29865748576503e-28*(pit_inv5*pit_inv)*(thetat_inv2*thetat_inv2) - 8.5671158651021397e-13*pit5*pit3*thetat_inv10 + 5.4200057337223301e-18*(pit5*pit3)*(thetat_inv10*thetat_inv2) + 8.58133791857099e-6*(pit10)*(thetat_inv5*thetat_inv3) + 2.6617045440598101e-14*(pit10)*(thetat_inv10*thetat_inv2) - 1.71242509570207e-37*thetat_inv10*thetat_inv2)
    V = main*main
    V *= V*0.0039
    return 1.0/V
    

def iapws97_region3_x(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.347826086956522e-08 - 0.91
    thetat = T*0.0015384615384615385 - 0.988
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    thetat2 = thetat*thetat
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    thetat_inv10 = thetat_inv5*thetat_inv5
    main = (-0.32085055136733398*pit*thetat_inv + 2.1578022250902002e-27*pit*(thetat_inv10) + thetat*(thetat*(thetat*(thetat*(thetat*(64866249228.068199*pit_inv + thetat*(-38264244845861000.0*pit2 + (thetat2)*((thetat2)*(1.6989448143359199e+21 + (thetat2)*((thetat2)*(-4.2599956229273801e+23*pit_inv2*pit_inv2 + 3.7737374129815101e+18*(pit_inv5*pit_inv3)) + 1.0731906585576701e+21*(pit_inv3)) - 1033632255988600.0*pit_inv5 - 5071008837229.1299*pit_inv5*pit_inv) - 3.2606864627931401e+20*pit3 + 2.6241320970635799e+24*(pit5*pit3)))) - 8515357334.8425798) + 39794900155318.398*(pit2*pit2)) - 0.00092472937839094499*pit_inv2*pit_inv2) + 1.8479081432077301e-6*(pit_inv2*pit_inv2) - 43235522531.974503*pit5 - 592874245598.60999*pit5*pit + 25818961427085.301*(pit5*pit3)) + 2.44200600688281 - 563199.25339166599*pit3 - 2.7538607767442098e-29*pit3*thetat_inv10*thetat_inv2 - 4.6230777187397299e-13*pit_inv3*thetat_inv2 + 16223.4569738433*(pit5)*(thetat_inv2) + 1.00824008584757e-7*(pit5)*(thetat_inv5*thetat_inv) + 1573381.97797544*(pit5*pit3)*(thetat_inv3) + 1.3306164728110601*(pit5*pit3)*(thetat_inv5*thetat_inv) - 0.092001193743114204*pit10*thetat_inv5*thetat_inv3 - 592910695.76253605*pit10*pit2*thetat_inv2*thetat_inv2 + 8470048.7061208691*(pit10*pit2)*(thetat_inv5) - 11.043375910954699*pit10*pit2*thetat_inv5*thetat_inv3 + 0.0022021376590542598*(pit10*pit2)*(thetat_inv10) + 4308676.5806146804*(pit10*pit2*pit2)*(thetat_inv5*thetat_inv) - 1192.28759669889*pit10*pit2*pit2*thetat_inv5*thetat_inv3 + 0.18133960351630199*(pit10*pit2*pit2)*(thetat_inv10) - 1.8302717326966e-5*pit10*pit2*pit2*thetat_inv10*thetat_inv2)
    V = main*0.0049
    return 1.0/V
    

def iapws97_region3_y(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.545454545454546e-08 - 0.996
    thetat = T*0.0015384615384615385 - 0.994
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit10 = pit5*pit5
    thetat2 = thetat*thetat
    thetat3 = thetat*thetat2
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    main = (thetat*(5834.4130522840696 + (thetat3)*(thetat*(thetat*((thetat2)*(-1.5909649090470799e+26*pit + 1.18973500934212e+25 - 2.66713136106469e+30*pit3) + 1.4933391705313e+27*(pit2*pit2)) - 13477896845792500.0 - 5.2711465785069603e+21*pit2) + 3.27777227273171e+18*(pit2) + 7.0510622439983402e+20*(pit3)) - 3818819062711000.0*pit5) + 496.21219715823901*thetat_inv*(pit2) - 3.15839902302021e-7*pit2*thetat_inv2*thetat_inv2 + 2.10017506281863e-17*(pit3)*(thetat_inv5*thetat_inv3) - 5.2559799502463301e-10*thetat_inv3 - 1.4537051255456199e-8*pit2*pit2*thetat_inv5*thetat_inv - 14979562.028764101*pit5*thetat_inv2 - 93780816955019.297*pit5*pit3*thetat_inv2 + 7.2466016558579696e-5*(pit5*pit3)*(thetat_inv5*thetat_inv3) + 5144114683.7638302*(pit10)*(thetat_inv5) - 82819.859404014103*pit10*pit2*thetat_inv5*thetat_inv3)
    V = main*main
    V *= V*0.0031
    return 1.0/V
    

def iapws97_region3_z(T, P):
    # This function was automatically generated. Do not edit it directly!
    pit = P*4.545454545454546e-08 - 0.993
    thetat = T*0.0015384615384615385 - 0.994
    pit_inv = 1.0/pit
    thetat_inv = 1.0/thetat
    pit2 = pit*pit
    pit3 = pit*pit2
    pit5 = pit3*pit2
    pit_inv2 = pit_inv*pit_inv
    pit_inv3 = pit_inv*pit_inv2
    pit_inv5 = pit_inv3*pit_inv2
    thetat2 = thetat*thetat
    thetat_inv2 = thetat_inv*thetat_inv
    thetat_inv3 = thetat_inv*thetat_inv2
    thetat_inv5 = thetat_inv3*thetat_inv2
    main = (-8.1973155961052001e-21*pit_inv*thetat_inv5*thetat_inv + thetat*(-62500479.117154002*pit + thetat*(thetat*(328380587890.65997 + (thetat2)*(thetat*(8.0319795746202006e+20*(pit2) + (thetat2)*(9238140070232500.0*(pit_inv2*pit_inv2) + 3277763028588600.0*(pit_inv5)) - 167170186672140.0*pit_inv3 - 3238999157299.6001*pit_inv2*pit_inv2 + 7288032747.7770996*(pit_inv5) - 4630574.3033124004*pit_inv5*pit_inv) + 663221436245.51001*(pit_inv3) - 1105981701.1840999*pit_inv2*pit_inv2) + 2.4400789229064998e-11*(pit_inv5*pit_inv3)) + 2537.4935870139002*(pit_inv2))) - 68285901137.457001*thetat_inv*pit5*pit - 3783.9104705594*pit3*thetat_inv2 - 2.0439701133835e-11*pit3*thetat_inv5*thetat_inv + 8.4225008041371002e-13*(pit_inv3)*(thetat_inv2) - 3739.6286292864002*pit5*pit*thetat_inv2*thetat_inv2 + 15.435572168146001*(pit5*pit)*(thetat_inv5) + 0.0097287654593861995*(pit5*pit)*(thetat_inv5*thetat_inv) + 3945360.4949706998*(pit5*pit3)*(thetat_inv2*thetat_inv2) - 0.00024848801561453999*pit5*pit3*thetat_inv5*thetat_inv3)
    V = main*main
    V *= V*0.0038
    return 1.0/V

def iapws97_region3_rho(T, P):
    r'''Calculate the mass density of water in region 3 of the IAPWS-97 standard.
    No cheking that the original point is in region 3 is performed.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]

    Returns
    -------
    rho : float
        Mass density of water in region 3, [kg/m^3]
        
    Notes
    -----
    Significant discontinuities exist between each region.

    Examples
    --------
    >>> iapws97_region3_rho(648.6, 22.5e6)
    353.06081088726
    '''
    region = iapws97_region_3(T, P)
    if region == REGION_3A:
        rho = iapws97_region3_a(T, P)
    elif region == REGION_3B:
        rho = iapws97_region3_b(T, P)
    elif region == REGION_3C:
        rho = iapws97_region3_c(T, P)
    elif region == REGION_3D:
        rho = iapws97_region3_d(T, P)
    elif region == REGION_3E:
        rho = iapws97_region3_e(T, P)
    elif region == REGION_3F:
        rho = iapws97_region3_f(T, P)
    elif region == REGION_3G:
        rho = iapws97_region3_g(T, P)
    elif region == REGION_3H:
        rho = iapws97_region3_h(T, P)
    elif region == REGION_3I:
        rho = iapws97_region3_i(T, P)
    elif region == REGION_3J:
        rho = iapws97_region3_j(T, P)
    elif region == REGION_3K:
        rho = iapws97_region3_k(T, P)
    elif region == REGION_3L:
        rho = iapws97_region3_l(T, P)
    elif region == REGION_3M:
        rho = iapws97_region3_m(T, P)
    elif region == REGION_3N:
        rho = iapws97_region3_n(T, P)
    elif region == REGION_3O:
        rho = iapws97_region3_o(T, P)
    elif region == REGION_3P:
        rho = iapws97_region3_p(T, P)
    elif region == REGION_3Q:
        rho = iapws97_region3_q(T, P)
    elif region == REGION_3R:
        rho = iapws97_region3_r(T, P)
    elif region == REGION_3S:
        rho = iapws97_region3_s(T, P)
    elif region == REGION_3T:
        rho = iapws97_region3_t(T, P)
    elif region == REGION_3U:
        rho = iapws97_region3_u(T, P)
    elif region == REGION_3V:
        rho = iapws97_region3_v(T, P)
    elif region == REGION_3W:
        rho = iapws97_region3_w(T, P)
    elif region == REGION_3X:
        rho = iapws97_region3_x(T, P)
    elif region == REGION_3Y:
        rho = iapws97_region3_y(T, P)
    elif region == REGION_3Z:
        rho = iapws97_region3_z(T, P)
    else:
        raise ValueError("Could not identify region 3 subregion")
    return rho


def iapws97_identify_region_TP(T, P):
    r'''Identify the main region given a temperature and pressure point
    according to the IAPWS097 standard.
    
    Raises a ValueError if the input point is out of bounds.

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]

    Returns
    -------
    region : int
        One of 1, 2, 3, or 5; region 4 is the two-phase region which cannot
        be reached with a TP call, [-]

    Examples
    --------
    >>> iapws97_identify_region_TP(400, 1e6)
    1
    '''
    under_623 = T <= 623.15
    if under_623:
        Psat = Psat_IAPWS(T)
    two_to_three = iapws97_boundary_2_3(T)
    if 273.15 <= T <= 623.15 and Psat < P <= 100E6:
        return 1
    elif 273.15 <= T <= 623.15 and P <= Psat:
        return 2
    elif 623.15 <= T <= 1073.15 and P <= 100E6 and P <= two_to_three:
        return 2
    elif 623.15 <= T <= 1073.15 and P <= 100E6 and P > two_to_three:
        return 3
    elif 1073.15 <= T <= 2273.15 and P <= 50E6:
        return 5
    else:
        raise ValueError("For box (1,2,3,4) 273.15 K <= T <= 1073.15 K and P <= 100 MPa;"
                         "for box 5, 1073.15 K <= T <= 2273.15 K and P <= 50 MPa.")

def iapws97_rho(T, P):
    r'''Calculate the mass density of water according to the IAPWS-97 standard.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]

    Returns
    -------
    rho : float
        Mass density of water, [kg/m^3]
        
    Notes
    -----
    Significant discontinuities exist between each region.

    Examples
    --------
    >>> iapws97_rho(648.6, 22.5e6)
    353.06081088726
    >>> iapws97_rho(330.0, 8e5)
    985.10498080770
    >>> iapws97_rho(823, 14e6)
    40.39293607288123
    >>> iapws97_rho(2000, 3e7)
    32.11456228328856
    '''
    R = 461.526
    region = iapws97_identify_region_TP(T, P)
    if region == 1:
        pi = P/16.53E6
        tau = 1386.0/T
        dG_dpi = iapws97_dG_dpi_region1(tau, pi)
        return P/(R*T*pi*dG_dpi)
    elif region == 2:
        pi = P/1E6
        tau = 540.0/T
        dG_dpi = 1.0/pi + iapws97_dGr_dpi_region2(tau, pi) 
        return P/(R*T*pi*dG_dpi)
    elif region == 3:
        return iapws97_region3_rho(T, P)
    elif region == 5:
        pi = P/1E6
        tau = 1000/T
        dG_dpi = 1.0/pi + iapws97_dGr_dpi_region5(tau, pi) 
        return P/(R*T*pi*dG_dpi)
    else:
        raise ValueError("Out of bounds")



def iapws95_dA_ddeltar(tau, delta):
    # Uses no constants from this file
    # Zeroth attempt Should try to replace each of the tau, delta powers by something else
    # Managed to make the derivative faster, so yeah it could have been better.
    tau05 = x0 = sqrt(tau)
    tau_quarter = sqrt(tau05)
    tau_eighth = sqrt(tau_quarter)
    x1 = tau
    x2 = delta
    delta2 = x3 = x23 = delta*delta
    delta3 = x4 = x31 = delta*x3 # delta 3
    delta4 = x16 = x49 = x4*delta #delta**4.0 # delta 4
    delta5 = x34 = delta*delta4
    delta6 = x18 = x36 = delta*delta5
    delta8 = x20 = delta6*delta2
    delta9 = x22 = delta*delta8
#    10, 11, 12, 13, 14 also used for delta but not needed
    x6 = delta - 1.0
    x7 = exp(-delta)
    x8 = x6*x7
    tau2 = x13 = tau*tau
    tau3 = x19 = tau*tau2
    tau4 = x5 = tau*tau3
    tau6 = x9 = tau3*tau3
    tau7 = tau6*tau
    tau9 = x17 = tau6*tau3
    tau10 = x29 = tau*tau9
    tau12 = tau10*tau2
    tau13 = x15 = tau10*tau3
    tau23 = tau10*tau13
    x10 = x1*x7
    x11 = x2*(delta - 2.0)
    x12 = x5*x7
    x14 = x4*x7*(delta - 4.0)
    x21 = x20*x5
    x24 = x23 + x23
    x25 = exp(-x23)
    x26 = tau*tau6*x25
    x27 = x2*(x23 - 1.0)
    x28 = x25*x27
    x30 = x25*x29
    x32 = exp(-x31)
    x33 = x4*(x23 - 2.0)
    x35 = x34*(x23 - 3.0)
    x37 = exp(-x36)
    x38 = x34*(x36 - 1.0)
    x39 = x37*x38
    x40 = tau23*tau23*tau4*x37
    x41 = x22*x25*(x23 - 5.0)
    x42 = tau4*tau4*x25
    x43 = 3.0*x31
    x44 = x4*(x43 - 4.0)
    x45 = tau23*x32
    x46 = x24 - 9.0
    x47 = x25*x46
    x48 = x20*x47
    x50 = delta - 1.0
    x51 = x50*x50
    dn1_2_23 = x51**(2.0/3.0)
    x52 = -20.0*x51
    x53 = x31*(-40.0*delta + 40.0 + 3.0/delta)
    c54 = (tau - 1.21)
    x54 = x53*exp(x52 - 150.0*c54*c54)
    x55 = -tau
    x51_2_x51sqrt = x51*x51*sqrt(x51)
    y50 =(0.32*x51*dn1_2_23 + x55 + 1.0)
    x56 = 0.2*x51*x51_2_x51sqrt + y50*y50
    x57 = 32.0*x51
    x58 = x61 = (tau - 1.0)*(tau - 1.0)
    x59 = 800.0*x58
    x60 = x51
    x62 = delta*x6
    x63 = 28.0*x51
    x64 = 700.0*x58
    x65 = x60**1.66666666666667
    c100 = (x55 + 0.32*x65 + 1.0)
    x66 = 0.2*x51_2_x51sqrt*x51 +  c100*c100
    c51 = exp(-x57 - x59)
    c52 = exp(-x63 - x64)
    x100 = (tau - 1.25)
    x101 = x66**(-0.05)
    x102 = x56**-0.05
    x67 = x62*(1.4*x51_2_x51sqrt + (-2.13333333333333*tau + 0.682666666666667*x65 + 2.13333333333333)*dn1_2_23)
    return (-3.6582165144204e-7*delta9*delta*x12*(delta - 11.0) + 3.277713668506e-5*delta9*delta2*x42*(x23 - 6)
            + 1.3251180074668e-12*delta9*delta3*x15*x7*(delta - 13.0) 
            + 0.00012537942082937*delta9*delta4*x29*(x49 + x49 - 7.0)*exp(-x49) 
            + 6.2639586912454e-10*delta9*delta5*x10*(delta - 15.0) - 0.0234599255063943*tau_quarter*tau_eighth*x3 
            - 0.52291067718716*tau_quarter*tau05*x2 + 7.8957634722828*tau_quarter*tau05*tau_eighth + 0.25709043003438*tau4*tau*x11*x7
            - 1.1537996422951e-9*tau*tau10*x22*x7*(delta - 10.0) + 6.6212605039687e-5*tau12*x8
            - 0.130840847171433*tau9*tau7*x3*x32*(x31 - 1.0) - 0.034994005463765*tau10*tau12*x32*x44 
            + 1.19434310126448*tau23*tau12*tau9*x39 - 1.90664983984428*tau23*tau23*x39 + 0.63605018690836*x0*x2
            - 0.035222982017504*x1*x28 + 0.0352357972408536*x1*x4 - 0.0083326504880713*x1*x48
            + 31.546140237781*x1*x54 - 8.7803203303561*x1 + 0.19232721156002*x10*x11 
            - 0.16074868486251*x12*x3*(delta - 3.0) + 0.040092828925807*x13*x14 + 0.029052336009585*x13*x48
            - 3.9343422603254e-7*x14*x15 + 7.5941377088144e-6*x16*x17*x7*(delta - 5.0)
            - 0.4780732991548*x16*x30*(x24 - 5.0) - 0.022446277332006*x16*x45*(x43 - 5.0)
            - 0.44264590335092*x17*x28 - 0.00031740616648314*x17*x41 - 0.00056250979351888*x18*x19*x7*(delta - 7.0)
            - 0.014180634400617*x18*x30*(x24 - 7.0) - 0.0099938293981612*x19*x25*x33 
            - 0.038615085574206*x19*x48 + 0.0016554050063734*x20*x42*x46 + 0.020393486513704*x21*x47
            + 1.5608652257135e-5*x21*x7*(delta - 9.0) - 0.041055881791896*x25*x35*x9 + 0.062717401425098*x26*x33 
            + 0.10793600908932*x26*(x24 - 1.0) + 0.80495339527056*x27*x30 - 0.58083399985759*x3*x30*(x24 - 3.0) 
            + 1.67133355696935e-9*x3*x40*(x36 + x36 - 1.0) + 1.48631859420682*x30*x33 + 0.27272870220686*x30*x35
            + 0.71047094555886*x38*x40 - 0.0039911143959082*x41*x9 + 0.076788197844621*x44*x45 
            - 2521.3154341695*x5*x53*exp(x52 - 250.0*x100*x100) + 0.66856572307965*x5*x8 
            - 31.306260323435*x54 - 0.14874640856724*x56*x102*x102*x102*(-56.0*x62*exp(-28.0*x60 - 700.0*x61)
            + c52) + 0.31806110878444*x56*x102*(-64.0*x62*exp(-32.0*x60 - 800.0*x61) 
            + c51) - 0.126434447282154*x101*x101*x101*x67*c52
            + 0.302158053345218*x101*x67*c51 - 0.20433810950965*x8*x9 + 0.012533547935523/x0)

            
            
def iapws95_d2A_d2deltar(tau, delta):
    # 4 sqrt; 4 exp; 2 power; 3 div; loads of multiplies and adds.
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta2*delta2
    x1 = delta - 1.0
    x2 = delta*x1
    x3 = x1*x1
    x4 = 800.0*delta2*x3 - 20.0*delta2 - 120.0*x2 + 3.0
    x5 = -20.0*x3
    x6 = (0.826446280991736*tau - 1.0)
    x6 = x4*exp(x5 - 219.615*x6*x6)
    x7 = delta*tau
    taurt = sqrt(tau)
    tau4rt = sqrt(taurt)
    tau8rt = sqrt(tau4rt)
    tau3 = tau*tau*tau
    tau4 = tau3*tau
    tau6 = tau3*tau3
    x9 = 2.0 - 2.0*delta
    x10 = (tau - 1.0)*(tau - 1.0)
    x11 = exp(-700.0*x10 - 28.0*x3)
    x3cbrt = x3**(1.0/3.0)
    x12 = x3*x3cbrt*x3cbrt
    x3rt = sqrt(x3)
    x16 = x3*x3*x3rt
    x13 = 0.32*x12 + 1.0 - tau
    x13 = 0.2*x3*x16 + x13*x13
    x14 = exp(-800.0*x10 - 32.0*x3)
    x20 = x13**(-0.05)
    x15 = x20*x20*x20
    x17 = x3cbrt*x3cbrt*(-2.1333333333333333*tau + 0.682666666666666755*x12 + 2.1333333333333333)
    x18 = 1.4*x16 + x17
    x19 = x1*x18
    x21 = (x16 + 0.714285714285714191*x17)
    x21 *= x3*x21
    x22 = (x18 + x3*((-2.84444444444444455*tau + 0.910222222222222488*x12 
                      + 2.84444444444444455)/x3cbrt + 2.275555555555556*x3*x3cbrt + 7.*x3*x3rt))
    x23 = exp(-delta)
    x25 = exp(-delta2)
    x58 = exp(-delta4)
    x59 = exp(-delta3)
    x24 = delta*x23
    x26 = delta2*x23
    delta5 = delta4*delta
    delta6 = delta4*delta2
    delta7 = delta*delta6
    delta8 = delta*delta7
    delta9 = delta*delta8
    delta12 = delta6*delta6
    delta15 = delta8*delta7
    x29 = delta2*x25
    x31 = x25*delta4
    x33 = x25*delta7
    x35 = x25*delta9
    x37 = x25*delta5*delta6
    x39 = x23*delta3
    x40 = x23*delta4
    x42 = x23*delta5
    x44 = x23*delta7
    x45 = x25*delta6
    x47 = x23*delta8
    x48 = x23*delta9
    x50 = x23*delta5*delta5
    x51 = x23*delta5*delta6
    x52 = x25*delta8
    x53 = x25*delta5*delta5
    x55 = x25*delta12
    x56 = delta*x25
    x57 = x25*delta3
    x60 = delta*x59
    x61 = delta4*x59
    x62 = delta7*x59
    x63 = (0.8*tau - 1.0)
    return (-0.046919851012788602*delta*tau8rt*tau4rt + 0.148746408567240002*delta*x11*(
            0.2499*x15/x13*x21 - 0.85*x15*x22) 
            - 0.318061108784439994*delta*x14*(0.0931*x20/x13*x21 
            - 0.949999999999999956*x20*x22) - 5042.63086833899979*delta*x4*tau4*exp(x5 - 390.625*x63*x63)
            - 62.6125206468700028*delta*x6 + taurt*(0.636050186908360016 - 0.52291067718716*tau4rt)
            + 0.105707391722560801*tau*delta2 - tau*(1.31543132516153401e-7*delta7*delta6*x23
            + tau*(-tau*(tau*(tau*(tau*(tau*(tau*(tau*(tau*(-0.00100303536663496002*delta15*delta5*x58 
            + 0.0077735240914209398*delta8*delta8*x58 + tau*(tau*(0.000132425210079373988*x23
            - 0.0000662126050396869941*x24 + x7*(tau3*(tau6*(tau*(tau3*tau6*tau6*tau6*(tau*tau*(11.439899039065681*delta15
            - 32.41304727735276*delta9 + 9.53324919922140168*delta3 + tau4*(-4.26282567335316021*delta15
            + 12.078006074500621*delta9 - 3.55235472779430017*delta3 + 3.67693382533256956e-8*delta6
            - 2.00560026836322003e-8*delta12 - 3.34266711393870005e-9)) - 7.16605860758687996*delta15
            + 20.3038327214961605*delta9 - 5.97171550632239967*delta3)*exp(-delta6) + 0.448925546640119966*delta2*x59
            - 0.80806598395221596*delta5*x59 + 0.20201649598805399*delta8*x59 - 0.921458374135452063*x60
            + 2.30364593533863005*x61 - 0.691093780601588992*x62) + 0.419928065565180031*x60 
            - 1.04982016391295008*x61 + 0.314946049173885023*x62) - 1.04672677737146391*delta3*x59
            + 0.392522541514298995*delta6*x59 + 0.261681694342865978*x59) - 1.32511800746680002e-12*x23*delta12
            + 4.72121071239048018e-6*x24 - 3.14747380826032012e-6*x26 + 3.93434226032540015e-7*x39
            - 2.06718409164820794e-10*x50 + 3.44530681941368034e-11*x51)) + 1.03841967806558997e-7*x47 
            - 2.30759928459020012e-8*x48 + 1.15379964229510002e-9*x50) - 7.59868993714932728*x25*delta5
            - 0.804953395270559979*x25 - 4.89314458888812087*x29 + 7.67603002421735958*x31 + 1.48687416460069*x33 
            + 0.0567225376024680025*x35 + 0.572835940275540079*x45 - 0.545457404413720037*x52
            - 0.0114095272954726698*delta12*x58 + 3.48500399914553993*x56 + 1.42978998508973909*x57)
            + 0.442645903350920022*x25 - 2.21322951675460011*x29 + 0.885291806701840045*x31
            - 0.00015188275417628801*x39 + 0.000075941377088144005*x40 - 7.59413770881439999e-6*x42 
            + 0.0142832774917412992*x52 - 0.0066655294961459402*x53 + 0.000634812332966279988*x55) 
            - 0.0000655542733701200008*x25*delta7*delta7 - 0.119189160458884807*x33 + 0.0629053902421892047*x35
            - 0.00662162002549360022*x37 - 0.00216329102121396019*x53 + 0.000819428417126499982*x55) 
            - 0.376304408550587999*x29 + 0.564456612825881998*x31 - 0.125434802850196009*x45 
            + 0.647616054535919972*x56 - 0.431744036357279981*x57) - 0.408676219019300013*x23 
            + 0.204338109509650007*x24 + 0.615838226878440032*x31 - 0.533726463294648013*x45 
            + 0.26171191139966099*x52 - 0.083813402314072194*x53 + 0.00798222879181639984*x55)
            - 0.51418086006875996*x23 + 1.02836172013751992*x24 - 0.25709043003437998*x26) 
            + 1.33713144615930002*x23 + 0.295926386095409999*x24 - 0.964492109175060008*x26 
            - 1.46833102898668777*x33 + 0.774952487520751965*x35 - 0.0815739460548159934*x37 
            + 0.160748684862510011*x39 - 0.00112382296251371978*x44 + 0.000280955740628429946*x47 
            + 0.0000246317294014894049*x48 - 8.04807633172487966e-6*x50 + 3.65821651442040008e-7*x51)
            - 0.00787513710926432062*x23*delta6 + 0.059962976388967193*x29 - 0.0899444645834507894*x31
            + 2.78028616134283224*x33 - 1.46737325181982814*x35 + 0.15446034229682401*x37 
            + 0.0236254113277929619*x42 + 0.000562509793518880044*x44 + 0.0199876587963223988*x45) 
            + 0.48111394710968397*x26 + 2.09176819269011993*x33 - 1.10398876836422999*x35 
            + 0.116209344038340004*x37 - 0.32074263140645598*x39 + 0.0400928289258069975*x40)
            - 1.87918760737362006e-8*x23*delta7*delta7 + 6.26395869124539993e-10*x23*delta15 + 0.384654423120040001*x23
            - 0.769308846240080002*x24 - 0.0352229820175039982*x25 + 0.192327211560020001*x26
            + 0.176114910087520005*x29 - 0.0704459640350079963*x31 - 0.599950835141133676*x33 
            + 0.316640718546709443*x35 - 0.0333306019522852034*x37) 
            - 8.32979887976543942*x11*x13*x15*(delta*(56.0*x3 - 1.0) + x9) 
            + 0.252868894564308*x11*x15*x19*(56.0*x2 - 1.0) 
            + 20.3559109622041596*x13*x20*x14*(delta*(64.0*x3 - 1.0) + x9) 
            - 0.604316106690436006*x14*x19*x20*(64.0*x2 - 1.0) + 63.0922804755619993*x6*x7)
