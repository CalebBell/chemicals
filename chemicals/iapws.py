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
from math import exp, log, sqrt, fsum
from chemicals.vapor_pressure import Psat_IAPWS, Tsat_IAPWS
from fluids.numerics import secant, newton, trunc_log, trunc_exp, horner, solve_2_direct, newton_system, damping_maintain_sign

__all__ = [
           'iapws97_boundary_2_3', 'iapws97_boundary_2_3_reverse',
           'iapws97_identify_region_TP', 'iapws97_region_3', 'iapws97_region3_rho',
           'iapws97_region1_rho', 'iapws97_region2_rho', 'iapws97_region5_rho',
           'iapws95_rho', 'iapws95_P',
           
           'iapws97_rho_extrapolated',
           'iapws97_rho', 'iapws97_P', 'iapws97_T', 
           
           'rhol_sat_IAPWS', 'rhog_sat_IAPWS', 'iapws95_Psat', 'rhol_sat_IAPWS95',
           'iapws95_saturation',
           
           'iapws95_A0', 'iapws95_dA0_dtau', 'iapws95_d2A0_dtau2', 'iapws95_d3A0_dtau3',
           
           'iapws95_Ar', 'iapws95_dAr_ddelta', 'iapws95_d2Ar_ddelta2', 'iapws95_d3Ar_ddelta3', 
           'iapws95_dAr_dtau', 'iapws95_d2Ar_dtau2',
           'iapws95_d2Ar_ddeltadtau',
           ]

__numba_additional_funcs__ = ['iapws97_region3_a', 'iapws97_region3_b', 'iapws97_region3_c', 
    'iapws97_region3_d', 'iapws97_region3_e', 'iapws97_region3_f', 'iapws97_region3_g', 'iapws97_region3_h',
    'iapws97_region3_i', 'iapws97_region3_j', 'iapws97_region3_k', 'iapws97_region3_l', 'iapws97_region3_m',
    'iapws97_region3_n', 'iapws97_region3_o', 'iapws97_region3_p', 'iapws97_region3_q', 'iapws97_region3_r', 
    'iapws97_region3_s', 'iapws97_region3_t', 'iapws97_region3_u', 'iapws97_region3_v', 'iapws97_region3_w', 
    'iapws97_region3_x', 'iapws97_region3_y', 'iapws97_region3_z',
    'iapws_97_Trho_err_region1', 'iapws_97_Trho_err_region2', 'iapws_97_Trho_err_region5',
    
    'iapws_97_Prho_err_region1', 'iapws_97_Prho_err_region2', 'iapws_97_Prho_err_region5',
    'iapws_97_Prho_err_region3',
    'iapws95_rho_err', 'iapws95_P_err',
    
    
   'iapws97_boundary_3uv', 'iapws97_boundary_3ef', 
   'iapws97_boundary_3ef', 'iapws97_boundary_3cd', 'iapws97_boundary_3gh',
   'iapws97_boundary_3ij', 'iapws97_boundary_3jk', 'iapws97_boundary_3mn', 
   'iapws97_boundary_3qu', 'iapws97_boundary_3rx', 'iapws97_boundary_3wx',
   'iapws97_boundary_3ab', 'iapws97_boundary_3op',

   
   'iapws97_G_region1', 'iapws97_dG_dpi_region1', 'iapws97_d2G_dpi2_region1',
   'iapws97_dG_dtau_region1', 'iapws97_d2G_d2tau_region1', 'iapws97_d2G_dpidtau_region1',
   
   'iapws97_Gr_region2', 'iapws97_dGr_dpi_region2', 'iapws97_d2Gr_d2pi_region2',
   'iapws97_dGr_dtau_region2', 'iapws97_d2Gr_d2tau_region2', 'iapws97_d2Gr_dpidtau_region2',
   'iapws97_G0_region2', 'iapws97_dG0_dtau_region2', 'iapws97_d2G0_d2tau_region2',
   
   'iapws97_A_region3', 'iapws97_dA_ddelta_region3', 'iapws97_d2A_d2delta_region3',
   'iapws97_dA_dtau_region3', 'iapws97_d2A_d2tau_region3', 'iapws97_d2A_ddeltadtau_region3',
   
   'iapws97_Gr_region5', 'iapws97_dGr_dpi_region5', 'iapws97_d2Gr_d2pi_region5', 
   'iapws97_dGr_dtau_region5', 'iapws97_d2Gr_d2tau_region5', 'iapws97_d2Gr_dpidtau_region5',
   'iapws97_G0_region5', 'iapws97_dG0_dtau_region5', 'iapws97_d2G0_d2tau_region5',
   

    ]

R95 = 461.51805 # Differs from the other formulation
R97 = 461.526

MW = 18.015268

Tc = 647.096
Tc_inv = 1.0/Tc

rhoc = 322.0
rhoc_inv = 1.0/rhoc

def use_mpmath_backend():
    import mpmath as mp
    globals()['exp'] = mp.exp
    globals()['log'] = mp.log
    globals()['sqrt'] = mp.sqrt
    globals()['R95'] = mp.mpf("461.51805")
    globals()['R97'] = mp.mpf("461.526")
    globals()['MW'] = mp.mpf("18.015268")
    globals()['Tc'] = mp.mpf("647.096")
    globals()['Tc_inv'] = 1/mp.mpf("647.096")
    globals()['rhoc'] = mp.mpf("322")
    globals()['rhoc_inv'] = 1/mp.mpf("322")
    
def reset_backend():
    import math
    globals()['exp'] = math.exp
    globals()['log'] = math.log
    globals()['sqrt'] = math.sqrt
    globals()['R95'] = 461.51805
    globals()['R97'] = 461.526
    globals()['MW'] = 18.015268
    globals()['Tc'] = 647.096
    globals()['Tc_inv'] = 1/647.096
    globals()['rhoc'] = 322.0
    globals()['rhoc_inv'] = 1/322.0


def iapws97_boundary_2_3(T):
    '''Above this pressure we are in region 3.
    
    >>> iapws97_boundary_2_3(0.623150000E3)
    16529164.2526216
    '''
    return (0.34805185628969E9 - T*(0.11671859879975E7 - 0.10192970039326E4*T))

def iapws97_boundary_2_3_reverse(P):
    '''
    >>> iapws97_boundary_2_3_reverse(16529164.2526216)
    623.15000000000
    '''
    return (116.85603437027114637*sqrt(7.1845068690062562659e-8*P - 1.0)
            + 572.54459862744727161)

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
#    T = 3.727888004*(P-22.064)+Tc
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

### Region 1
    
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
    tmp = pit3*pit3*pit3
    pit18 = tmp*tmp
    pit27 = pit18*tmp
    
    taut3 = taut*taut*taut
    taut7 = taut3*taut3*taut

    taut_inv2 = taut_inv*taut_inv
    taut_inv3 = taut_inv*taut_inv2
    tmp = taut_inv3*taut_inv3 # 6
    tmp *= tmp  # 12; also temporary variable
    taut_inv29 = tmp*tmp*taut_inv2 # 26
    taut_inv38 = taut_inv29*tmp # 26 + 12 = 38 end
    taut_inv29 *= taut_inv3 # 29 end    
    # some savings remain below via horner with appropriate testing
    return (-0.00001696247879118719984*pit - 0.0001900778690703240094*pit*taut_inv3*taut_inv
            - 5.112307687206180469e-9*pit*taut3*taut3 + 0.00009532278781397400275*taut
            - 0.0006000356158605199813 - 7.820546747472120262e-6*pit2*taut_inv2
            - 0.00002691033828959999829*pit2*taut_inv3*taut_inv2 - 1.721007592550879906e-12*pit2*taut7*taut3 
            - 8.103399372023399331e-6*pit3*taut_inv3*taut_inv3*taut_inv2 - 8.828369066169199338e-6*taut3
            - 0.0009436864214653399542*taut_inv3 - 9.757927889155040732e-9*pit3*pit3*taut_inv3*taut_inv3 
            - 7.131208975318960007e-8*pit3*pit3*taut_inv3*taut_inv3*taut_inv3*taut_inv2 
            - 1.453899925951880029e-15*taut7*taut7*taut3 - 2.888009514412301439e-16*pit18*pit*taut_inv29 
            + 7.326023761231625652e-18*pit18*pit3*taut_inv29*taut_inv2 + 2.138465471018954057e-20*pit27*taut_inv38
            - 1.039443169686176943e-20*pit27*pit*taut_inv38*taut_inv
            + 1.695212796070572203e-21*pit27*pit2*taut_inv38*taut_inv2 
            - 9.278879059411833807e-23*pit27*pit3*taut_inv38*taut_inv3)


def iapws97_dG_dtau_region1(tau, pi):
    pit = 7.1 - pi
    taut = tau - 1.222
    taut_inv = 1.0/taut
    
    pit2 = pit*pit
    pit4 = pit2*pit2
    tmp = pit4*pit4 # 8
    pit21 = tmp*tmp*pit4*pit # 20
    pit29 = pit21*tmp
    
    taut2 = taut*taut
    taut4 = taut2*taut2
    
    taut_inv2 = taut_inv*taut_inv
    taut_inv3 = taut_inv2*taut_inv
    tmp = taut_inv3*taut_inv3 # 6
    taut_inv9 = tmp*taut_inv3
    taut_inv30 = taut_inv9*tmp # 15
    taut_inv30 *= taut_inv30
    return (-0.02184171717541399937*pit - 0.0001585150739097900138*pit*taut2 +
            0.01899006821841900047*pit*taut_inv2 + 0.004249441109611180011*pit*taut_inv3*taut_inv3*taut_inv2
            - 0.002548717211142359929*pit*taut_inv9*taut_inv - 1.915839267757440068*taut
            + 3.385516916838500201 + 0.00004766139390698700138*pit2 - 0.00001324255359925379985*pit2*taut2 
            + 0.001415529632198009877*pit2*taut_inv3*taut_inv 
            - 1.235814937059098064e-14*pit2*taut4*taut4*taut4*taut4 + 0.4731611553968400052*taut2 
            + 0.8454818716911399745*taut_inv2 - 5.112307687206180469e-9*pit2*pit*taut4*taut
            + 0.0001267185793802160063*pit2*pit*taut_inv3*taut_inv2 - 0.06646566879800400174*taut2*taut
            - 0.2926594242633400178*taut_inv3 + 1.303424457912020044e-6*pit4*taut_inv3
            + 0.00001121264095399999929*pit4*taut_inv3*taut_inv3 - 1.434172993792399922e-12*pit4*taut4*taut4*taut 
            + 0.00406073149917839972*taut4 + 3.241359748809359987e-6*pit4*pit*taut_inv9 
            + 1.045492273838040137e-9*pit4*pit4*taut_inv3*taut_inv3*taut_inv 
            + 1.400773191580509978e-8*pit4*pit4*taut_inv9*taut_inv3 + 1.99410180757039883e-17*pit21*taut_inv30 
            - 4.488275426841510012e-19*pit21*pit2*taut_inv30*taut_inv2 
            - 1.000759703186210033e-21*pit29*taut_inv30*taut_inv9 
            + 4.659572829627690075e-22*pit29*pit*taut_inv30*taut_inv9*taut_inv 
            - 7.291237832561599838e-23*pit29*pit2*taut_inv30*taut_inv9*taut_inv2 
            + 3.835020578990777884e-24*pit29*pit2*pit*taut_inv30*taut_inv9*taut_inv3)

def iapws97_d2G_d2tau_region1(tau, pi):
    pit = 7.1 - pi
    taut = tau - 1.222
    taut_inv = 1.0/taut
    
    pit2 = pit*pit
    pit4 = pit2*pit2
    taut_inv2 = pit4*pit4 # abuse taut_inv2 variable as a temporary
    pit21 = taut_inv2*taut_inv2*pit4*pit
    pit29 = pit21*taut_inv2
    
    taut2 = taut*taut
    taut4 = taut2*taut2
    
    taut_inv2 = taut_inv*taut_inv
    taut_inv4 = taut_inv2*taut_inv2
    taut_inv9 = taut_inv4*taut_inv4*taut_inv
    taut_inv31 = taut_inv9*taut_inv9*taut_inv9*taut_inv4
    return (-0.0003170301478195800275*pit*taut - 0.03798013643683800095*pit*taut_inv2*taut_inv
            - 0.03399552887688944008*pit*taut_inv9 + 0.02548717211142359843*pit*taut_inv9*taut_inv2
            + 0.9463223107936800105*taut - 0.00002648510719850759971*taut*pit2 - 1.915839267757440068
            - 0.005662118528792039508*pit2*taut_inv4*taut_inv
            - 1.977303899294556903e-13*pit2*taut4*taut4*taut4*taut2*taut - 0.1993970063940120052*taut2
            - 2.556153843603090152e-8*pit2*pit*taut4 - 0.0006335928969010799772*pit2*pit*taut_inv4*taut_inv2
            + 0.01624292599671359888*taut2*taut - 1.690963743382279949*taut_inv2*taut_inv 
            - 3.910273373736060131e-6*pit4*taut_inv4 - 0.00006727584572399999572*pit4*taut_inv4*taut_inv2*taut_inv 
            - 1.29075569441316001e-11*pit4*taut4*taut4 + 0.8779782727900200534*taut_inv4 
            - 0.0000291722377392842403*pit4*pit*taut_inv9*taut_inv - 7.318445916866280962e-9*pit4*pit4*taut_inv4*taut_inv4 
            - 1.680927829896611973e-7*pit4*pit4*taut_inv9*taut_inv4 - 5.982305422711196613e-16*pit21*taut_inv31 
            + 1.436248136589283204e-17*pit21*pit2*taut_inv31*taut_inv2 + 3.902962842426219073e-20*pit29*taut_inv31*taut_inv9
            - 1.863829131851076105e-20*pit29*pit*taut_inv31*taut_inv9*taut_inv
            + 2.989407511350256051e-21*pit29*pit2*taut_inv31*taut_inv9*taut_inv2 
            - 1.610708643176126667e-22*pit29*pit2*pit*taut_inv31*taut_inv9*taut_inv2*taut_inv)

def iapws97_d2G_dpidtau_region1(tau, pi):
    pit = 7.1 - pi
    taut = tau - 1.222
    taut_inv = 1.0/taut
    
    pit2 = pit*pit
    pit3 = pit2*pit
    pit7 = pit3*pit3*pit
    pit28 = pit7*pit7
    pit28 *= pit28
    
    taut2 = taut*taut
    taut8 = taut2*taut2
    taut8 *= taut8
    
    taut_inv2 = taut_inv*taut_inv
    taut_inv3 = taut_inv2*taut_inv
    taut_inv30 = taut_inv3*taut_inv3 # 6 temporarily
    taut_inv9 = taut_inv30*taut_inv3
    taut_inv30 = taut_inv9*taut_inv30
    taut_inv30 *= taut_inv30 # 30 finally
    return (-0.00009532278781397400275*pit + 0.00002648510719850759971*pit*taut2 - 0.002831059264396019754*pit*taut_inv3*taut_inv + 2.471629874118196128e-14*pit*taut8*taut8 + 0.02184171717541399937 + 1.533692306161854223e-8*pit2*taut2*taut2*taut - 0.0003801557381406480188*pit2*taut_inv3*taut_inv2 + 0.0001585150739097900138*taut2 - 0.01899006821841900047*taut_inv2 - 5.213697831648080174e-6*pit3*taut_inv3 - 0.00004485056381599999715*pit3*taut_inv3*taut_inv3 + 5.736691975169599686e-12*pit3*taut8*taut - 0.00001620679874404679866*pit3*pit*taut_inv9 - 8.3639381907043211e-9*pit7*taut_inv3*taut_inv3*taut_inv - 1.120618553264407982e-7*pit7*taut_inv9*taut_inv3 - 0.004249441109611180011*taut_inv3*taut_inv3*taut_inv2 + 0.002548717211142359929*taut_inv9*taut_inv - 4.187613795897837728e-16*pit7*pit7*pit3*pit3*taut_inv30 + 1.03230334817354737e-17*pit7*pit7*pit7*pit*taut_inv30*taut_inv2 + 2.902203139240009378e-20*pit28*taut_inv30*taut_inv9 - 1.397871848888307079e-20*pit28*pit*taut_inv30*taut_inv9*taut_inv + 2.26028372809409602e-21*pit28*pit2*taut_inv30*taut_inv9*taut_inv2 - 1.227206585277048923e-22*pit28*pit3*taut_inv30*taut_inv9*taut_inv3)

### Region 2

def iapws97_G0_region2(tau, pi):
    tau_inv = 1.0/tau
    return (tau*(tau*(0.0212684637533070015*tau - 0.284086324607719987) 
            + 10.0866559680180004) + tau_inv*(tau_inv*(tau_inv*(tau_inv*(0.0714527380814549973 
            - 0.00560879112830200022*tau_inv) - 0.407104982239279989) + 1.42408191714439991)
            - 4.38395113194500041) + log(pi) - 9.69276865002169963)

def iapws97_dG0_dtau_region2(tau, pi):
    # does not depend on pi but leave as argument for consistency
    tau_inv = 1.0/tau
    return (tau*(0.0638053912599210044*tau - 0.568172649215439973)
            + tau_inv*tau_inv*(tau_inv*(tau_inv*(tau_inv*(0.0280439556415100003*tau_inv
            - 0.285810952325819989) + 1.22131494671783991) - 2.84816383428879982) 
            + 4.38395113194500041) + 10.0866559680180004)

def iapws97_d2G0_d2tau_region2(tau, pi):
    # does not depend on pi but leave as argument for consistency
    tau_inv = 1.0/tau
    return (0.127610782519842009*tau + tau_inv*tau_inv*tau_inv*(tau_inv*(tau_inv*(tau_inv*(1.42905476162910006 
            - 0.168263733849060015*tau_inv) - 4.88525978687135964) + 8.54449150286639991)
            - 8.76790226389000082) - 0.568172649215439973)



def iapws97_Gr_region2(tau, pi):
    taut = tau - 0.5
    pi2 = pi*pi
    pi3 = pi2*pi
    pi4 = pi2*pi2
    pi8 = pi4*pi4
    pi20 = pi8*pi8*pi4
    
    taut2 = taut*taut
    taut3 = taut2*taut
    taut4 = taut2*taut2
    taut10 = taut4*taut4*taut2
    taut25 = taut10*taut10*taut2*taut3
    return (-0.01783486229235799886*pi*taut - 0.001773174247321299916*pi - 0.04599601369636500264*pi*taut2 - 0.05758125908343200011*pi*taut3 - 0.05032527872793000207*pi*taut4*taut2 - 0.00003303264167020299999*taut*pi2 + 4.387066728443500103e-7*taut*pi3 - 7.884730955936700091e-10*taut*pi4 - 0.000189489875163149999*pi2*taut2 - 0.003939277724335500143*pi2*taut4 - 0.04379729565057299823*pi2*taut4*taut3 - 0.00002667454791408700141*pi2*taut25*taut10*taut + 2.048173769230899887e-8*pi3 - 0.0000322776772385700023*pi3*taut3 - 0.001503392454214799983*pi3*taut4*taut2 - 0.04066825356264899827*pi3*taut25*taut10 + 1.279071785228500082e-8*pi4*taut2 + 4.82253727185070016e-7*pi4*taut3 + 2.292207633766100113e-6*pi4*pi*taut4*taut3 - 1.671476645106100115e-11*pi4*pi2*taut3 - 0.002117147232135499837*pi4*pi2*taut10*taut4*taut2 - 23.89574193410399872*pi4*pi2*taut25*taut10 - 5.905956432427000368e-18*pi4*pi3 - 1.262180889910100042e-6*pi4*pi3*taut10*taut - 0.03894684243573900279*pi4*pi3*taut25 + 1.125621136045899983e-11*pi8*taut4*taut4 - 8.231134089799800435*pi8*taut25*taut10*taut + 1.980971280208800021e-8*pi8*pi*taut10*taut3 + 1.040696521017399955e-19*pi8*pi2*taut4 - 1.023474709592900015e-13*pi8*pi2*taut10 - 1.001817937951099974e-9*pi8*pi2*taut10*taut4 - 8.088290864698499771e-11*pi8*pi8*taut25*taut4 + 0.1069303187940899985*pi8*pi8*taut25*taut25 - 0.3366225057417099875*pi8*pi8*pi2*taut25*taut25*taut4*taut3 + 8.918584535542099871e-25*pi20*taut10*taut10 + 3.062931687623199748e-13*pi20*taut25*taut10 - 4.200246769820800092e-6*pi20*taut25*taut10*taut10*taut3 - 5.905602968563900274e-26*pi20*pi*taut10*taut10*taut + 3.782694761345700151e-6*pi20*pi2*taut25*taut25*taut3 - 1.276860893468100005e-15*pi20*pi3*taut25*taut10*taut4 + 7.308761059506100019e-29*pi20*pi4*taut25*taut + 5.541471535077800115e-17*pi20*pi4*taut25*taut10*taut4*taut - 9.436970724120999844e-7*pi20*pi4*taut25*taut25*taut4*taut4)

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
    >>> iapws97_dGr_dpi_region2(.656, 16.0)
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

def iapws97_d2Gr_d2pi_region2(tau, pi):
    taut = tau - 0.5
    pi2 = pi*pi
    pi4 = pi2*pi2
    pi8 = pi4*pi4
    pi18 = pi2*pi8*pi8
    
    taut2 = taut*taut
    taut3 = taut2*taut
    taut4 = taut2*taut2
    taut10 = taut4*taut4*taut2
    taut25 = taut10*taut10*taut4*taut
    return (2.632240037066100168e-6*pi*taut + 1.228904261538539999e-7*pi - 0.0001936660634314200003*pi*taut3 - 0.00902035472528879903*pi*taut4*taut2 - 0.2440095213758939896*pi*taut25*taut10 - 0.00006606528334040599997*taut - 9.461677147124039696e-9*taut*pi2 + 1.534886142274200165e-7*pi2*taut2 + 5.787044726220840615e-6*pi2*taut3 - 0.0003789797503262999981*taut2 + 0.00004584415267532200057*pi2*pi*taut4*taut3 - 5.014429935318299763e-10*pi4*taut3 - 0.06351441696406499859*pi4*taut10*taut4*taut2 - 716.8722580231200254*pi4*taut25*taut10 - 0.007878555448671000286*taut4 - 2.480501701619340401e-16*pi4*pi - 0.00005301159737622419582*pi4*pi*taut10*taut - 1.635767382301038353*pi4*pi*taut25 + 6.303478361857040293e-10*pi4*pi2*taut4*taut4 - 460.9435090287888102*pi4*pi2*taut25*taut10*taut + 1.426299321750336068e-6*pi4*pi2*pi*taut10*taut3 - 0.08759459130114599645*taut4*taut3 + 9.366268689156599206e-18*pi8*taut4 - 9.211272386336099881e-12*pi8*taut10 - 9.016361441559899391e-8*pi8*taut10*taut4 - 1.941189807527639966e-8*pi8*pi4*pi2*taut25*taut4 + 25.66327651058159987*pi8*pi4*pi2*taut25*taut25 - 103.0064867569632554*pi8*pi8*taut25*taut25*taut4*taut3 + 3.389062123505998002e-22*pi18*taut10*taut10 + 1.16391404129681584e-10*pi18*taut25*taut10 - 0.001596093772531903933*pi18*taut25*taut10*taut10*taut3 - 2.4803532467968384e-23*pi18*pi*taut10*taut10*taut + 0.001747604979741713581*pi18*pi2*taut25*taut25*taut3 - 6.460916120948586575e-13*pi18*pi2*pi*taut25*taut10*taut4 + 4.034436104847367623e-26*pi18*pi4*taut25*taut + 3.058892287362945313e-14*pi18*pi4*taut25*taut10*taut4*taut - 0.0005209207839714791177*pi18*pi4*taut25*taut25*taut4*taut4 - 0.00005334909582817400282*taut25*taut10*taut)

def iapws97_dGr_dtau_region2(tau, pi):
    taut = tau - 0.5
    pi2 = pi*pi
    pi3 = pi2*pi
    pi4 = pi2*pi2
    pi8 = pi4*pi4
    pi20 = pi8*pi8*pi4
    
    taut2 = taut*taut
    taut4 = taut2*taut2
    taut5 = taut4*taut
    taut13 = taut4*taut4*taut5 # 13
    taut34 = taut13*taut4 # 9 + 8 = 17
    taut34 *= taut34 # 34 end
    return (-0.09199202739273000529*pi*taut - 0.01783486229235799886*pi - 0.1727437772502959934*pi*taut2 - 0.3019516723675800263*pi*taut5 - 0.0003789797503262999981*taut*pi2 + 2.558143570457000165e-8*taut*pi4 - 0.00003303264167020299999*pi2 - 0.01575711089734200057*pi2*taut2*taut - 0.3065810695540109876*pi2*taut5*taut - 0.0009602837249071320778*pi2*taut34*taut + 4.387066728443500103e-7*pi3 - 0.00009683303171571000013*pi3*taut2 - 0.00902035472528879903*pi3*taut5 - 1.423388874692715023*pi3*taut34 - 7.884730955936700091e-10*pi4 + 1.446761181555210154e-6*pi4*taut2 + 0.00001604545343636269952*pi4*pi*taut5*taut - 5.014429935318300022e-11*pi4*pi2*taut2 - 0.0338743557141679974*pi4*pi2*taut13*taut2 - 836.3509676936399728*pi4*pi2*taut34 - 0.00001388398978901110004*pi4*pi3*taut5*taut5 - 0.9736710608934751043*pi4*pi3*taut13*taut5*taut5*taut + 9.004969088367199864e-11*pi8*taut5*taut2 - 296.3208272327927943*pi8*taut34*taut + 2.575262664271440094e-7*pi8*pi*taut5*taut5*taut2 + 4.162786084069599818e-19*pi8*pi2*taut2*taut - 1.023474709592899964e-12*pi8*pi2*taut5*taut2*taut2 - 1.402545113131540004e-8*pi8*pi2*taut13 - 2.345604350762564972e-9*pi8*pi8*taut13*taut13*taut2 + 5.34651593970450012*pi8*pi8*taut34*taut13*taut2 - 19.18748282727747068*pi8*pi8*pi2*taut34*taut13*taut5*taut2*taut2 + 1.783716907108420048e-23*pi20*taut13*taut5*taut + 1.072026090668119912e-11*pi20*taut34 - 0.0002016118449513984044*pi20*taut34*taut13 - 1.240176623398419126e-24*pi20*pi*taut13*taut5*taut2 + 0.0002004828223513221118*pi20*pi2*taut34*taut13*taut5 - 4.979757484525589725e-14*pi20*pi3*taut34*taut2*taut2 + 1.90027787547158596e-27*pi20*pi4*taut13*taut5*taut5*taut2 + 2.216588614031120095e-15*pi20*pi4*taut34*taut5 - 0.00005473443019990179931*pi20*pi4*taut34*taut13*taut5*taut5)

def iapws97_d2Gr_d2tau_region2(tau, pi):
    taut = tau - 0.5
    pi2 = pi*pi
    pi3 = pi2*pi
    pi4 = pi2*pi2
    pi8 = pi4*pi4
    pi20 = pi8*pi8*pi4
    
    taut2 = taut*taut
    taut4 = taut2*taut2
    taut5 = taut4*taut
    taut18 = taut4*taut4
    taut18 *= taut18 # 16
    taut33 = taut18*taut18*taut #32
    taut18 *= taut2 # 18
    return (-0.3454875545005919868*pi*taut - 0.09199202739273000529*pi - 1.509758361837900242*pi*taut4 - 0.0001936660634314200003*taut*pi3 + 2.893522363110420308e-6*taut*pi4 - 1.002885987063660004e-10*taut*pi4*pi2 - 0.0003789797503262999981*pi2 - 0.04727133269202600518*pi2*taut2 - 1.839486417324065926*pi2*taut5 - 0.03360993037174962034*pi2*taut33*taut - 0.04510177362644399168*pi3*taut4 - 48.3952217395523121*pi3*taut33 + 2.558143570457000165e-8*pi4 + 0.00009627272061817619036*pi4*pi*taut5 - 0.5081153357125199888*pi4*pi2*taut5*taut5*taut4 - 28435.93290158375748*pi4*pi2*taut33 - 0.0001388398978901110071*pi4*pi3*taut5*taut4 - 23.3681054614434025*pi4*pi3*taut18*taut5 + 6.303478361857040293e-10*pi8*taut5*taut - 10371.22895314774723*pi8*taut33*taut + 3.090315197125728325e-6*pi8*pi*taut5*taut5*taut + 1.248835825220879946e-18*pi8*pi2*taut2 - 9.211272386336099881e-12*pi8*pi2*taut5*taut2*taut - 1.823308647071001956e-7*pi8*pi2*taut5*taut5*taut2 - 6.567692182135182584e-8*pi8*pi8*taut18*taut5*taut4 + 261.9792810455205085*pi8*pi8*taut33*taut5*taut5*taut5 - 1074.499038327538301*pi8*pi8*pi2*taut33*taut18*taut4 + 3.389062123505998002e-22*pi20*taut18 + 3.644888708271607926e-10*pi20*taut33 - 0.009475756712715725089*pi20*taut33*taut5*taut5*taut2*taut - 2.4803532467968384e-23*pi20*pi*taut18*taut + 0.01042510676226874981*pi20*pi2*taut33*taut18 - 1.892307844119724095e-12*pi20*pi3*taut33*taut4 + 4.750694688678964685e-26*pi20*pi4*taut18*taut5*taut + 8.644695594721368726e-14*pi20*pi4*taut33*taut5 - 0.003119862521394402635*pi20*pi4*taut33*taut18*taut5)

def iapws97_d2Gr_dpidtau_region2(tau, pi):
    taut = tau - 0.5
    pi2 = pi*pi
    pi3 = pi2*pi
    pi6 = pi3*pi3
    pi17 = pi6*pi6*pi3*pi2

    taut2 = taut*taut
    taut4 = taut2*taut2
    taut5 = taut4*taut
    taut13 = taut4*taut4*taut5 # 13
    taut34 = taut13*taut4 # 9 + 8 = 17
    taut34 *= taut34 # 34 end
    return (-0.0007579595006525999962*pi*taut - 0.00006606528334040599997*pi - 0.03151422179468400114*pi*taut2*taut - 0.6131621391080219752*pi*taut5*taut - 0.001920567449814264156*pi*taut34*taut - 0.09199202739273000529*taut + 1.023257428182800066e-7*taut*pi3 - 0.01783486229235799886 + 1.316120018533050084e-6*pi2 - 0.0002904990951471300275*pi2*taut2 - 0.02706106417586639709*pi2*taut5 - 4.270166624078144402*pi2*taut34 - 0.1727437772502959934*taut2 - 3.153892382374680036e-9*pi3 + 5.787044726220840615e-6*pi3*taut2 + 0.000080227267181813501*pi3*pi*taut5*taut - 3.008657961190980272e-10*pi3*pi2*taut2 - 0.2032461342850079844*pi3*pi2*taut13*taut2 - 5018.10580616183961*pi3*pi2*taut34 - 0.3019516723675800263*taut5 - 0.00009718792852307769686*pi6*taut5*taut5 - 6.815697426254326174*pi6*taut13*taut5*taut5*taut + 7.203975270693759891e-10*pi6*pi*taut5*taut2 - 2370.566617862342355*pi6*pi*taut34*taut + 2.317736397844296243e-6*pi6*pi2*taut5*taut5*taut2 + 4.162786084069599818e-18*pi6*pi3*taut2*taut - 1.023474709592900005e-11*pi6*pi3*taut5*taut2*taut2 - 1.402545113131539905e-7*pi6*pi3*taut13 - 3.752966961220103956e-8*pi6*pi6*pi3*taut13*taut13*taut2 + 85.54425503527200192*pi6*pi6*pi3*taut34*taut13*taut2 - 345.3746908909944295*pi17*taut34*taut13*taut5*taut2*taut2 + 3.567433814216839978e-22*pi17*pi2*taut13*taut5*taut + 2.144052181336239759e-10*pi17*pi2*taut34 - 0.004032236899027968197*pi17*pi2*taut34*taut13 - 2.604370909136680202e-23*pi17*pi3*taut13*taut5*taut2 + 0.004410622091729086459*pi17*pi3*pi*taut34*taut13*taut5 - 1.145344221440885637e-12*pi17*pi3*pi2*taut34*taut2*taut2 + 4.560666901131806878e-26*pi17*pi6*taut13*taut5*taut5*taut2 + 5.319812673674687913e-14*pi17*pi6*taut34*taut5 - 0.001313626324797643238*pi17*pi6*taut34*taut13*taut5*taut5)


### Region 3 A formulations
def iapws97_A_region3(tau, delta):
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta3*delta
    delta5 = delta4*delta
    
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau6 = tau4*tau2
    tau13 = tau6*tau6*tau
    return (-1.265431547771399989*delta*tau2 - 1.152440780668100073*delta*tau6 + 0.8852104398431800414*delta*tau13*tau2 - 0.6420776518160700164*delta*tau13*tau4 + 20.94439697430700065*tau + 0.1077051262633200029*tau*delta5 - 0.0001655767979503699929*tau*delta5*delta5 + 1.065807002851300034*log(delta) - 15.73284529023900014 + 0.3849346018667100244*delta2 - 0.8521470882420599802*delta2*tau2 + 4.897228154187700078*delta2*tau6 - 3.050261725696500115*delta2*tau6*tau + 0.03942053687915399868*delta2*tau13*tau6*tau2*tau + 0.1255840842430800131*delta2*tau13*tau13 - 7.686770787871600064*tau2 - 0.2799932969871000155*delta3 + 1.389979956946000073*delta3*tau2 - 2.018991502357000201*delta3*tau4 - 0.008214763717396300277*delta3*tau13*tau2*tau - 0.4759603573492299788*delta3*tau13*tau13 + 0.0439840744735000011*delta4 - 0.4447643542873899736*delta4*tau2 + 0.9057207071973299994*delta4*tau4 + 0.7052245008796700354*delta4*tau13*tau13 - 0.3291362325895400009*delta5*tau2*tau - 0.5087106204115799946*delta5*tau13*tau13 - 0.02217540087309599964*delta5*delta + 0.0942607516650919991*delta5*delta*tau2 + 0.1643627844796100024*delta5*delta*tau13*tau13 - 0.01350337224134800021*delta5*delta2*tau2 + 2.618594778795400035*tau6*tau - 0.01483434535247200002*delta5*delta3*tau13*tau13 + 0.0005792295362808399465*delta5*delta4*tau2 + 0.00323089047037109986*delta5*delta4*tau13*tau13 + 0.00008096480299621500545*delta5*delta5 - 2.808078114861999985*tau6*tau4 - 0.00004492389906181499664*delta5*delta5*delta*tau13*tau13 + 1.205336969651700008*tau6*tau6 - 0.008456681281250200827*tau13*tau6*tau4)
    
def iapws97_dA_ddelta_region3(tau, delta):
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta3*delta
    delta5 = delta4*delta
    
    tau2 = tau*tau
    tau3 = tau2*tau
    tau6 = tau3*tau3
    tau13 = tau6*tau6*tau
    return (0.7698692037334200489*delta - 1.70429417648411996*delta*tau2 + 9.794456308375400155*delta*tau6 - 6.100523451393000229*delta*tau6*tau + 0.07884107375830799735*delta*tau13*tau6*tau3 + 0.2511681684861600261*delta*tau13*tau13 + 0.5385256313166000286*tau*delta4 - 0.001655767979503699984*tau*delta5*delta4 + 1.065807002851300034/delta - 0.8399798909613001019*delta2 + 4.169939870838000218*delta2*tau2 - 6.056974507071000602*delta2*tau3*tau - 0.0246442911521888991*delta2*tau13*tau3 - 1.427881072047689992*delta2*tau13*tau13 - 1.265431547771399989*tau2 + 0.1759362978940000044*delta3 - 1.779057417149559894*delta3*tau2 + 3.622882828789319998*delta3*tau3*tau + 2.820898003518680142*delta3*tau13*tau13 - 1.645681162947699949*delta4*tau3 - 2.543553102057900084*delta4*tau13*tau13 - 0.1330524052385760048*delta5 + 0.5655645099905519668*delta5*tau2 + 0.9861767068776600142*delta5*tau13*tau13 - 0.09452360568943600494*delta5*delta*tau2 - 1.152440780668100073*tau6 - 0.1186747628197760002*delta5*delta2*tau13*tau13 + 0.005213065826527559302*delta5*delta3*tau2 + 0.02907801423333989874*delta5*delta3*tau13*tau13 + 0.0008096480299621500003*delta5*delta4 - 0.0004941628896799649291*delta5*delta5*tau13*tau13 + 0.8852104398431800414*tau13*tau2 - 0.6420776518160700164*tau13*tau3*tau)
    
def iapws97_d2A_d2delta_region3(tau, delta):
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta3*delta
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau6 = tau4*tau2
    tau20 = tau6*tau6*tau6*tau2
    return ( - 1.065807002851300034/delta2 -1.679959781922600204*delta + 8.339879741676000435*delta*tau2 - 12.1139490141420012*delta*tau4 - 0.04928858230437779819*delta*tau6*tau6*tau4 - 2.855762144095379984*delta*tau20*tau6 + 2.154102525266400114*tau*delta3 - 0.0149019118155332992*tau*delta4*delta4 + 0.7698692037334200489 + 0.5278088936820000132*delta2 - 5.337172251448679461*delta2*tau2 + 10.86864848636795955*delta2*tau4 + 8.462694010556040425*delta2*tau20*tau6 - 1.70429417648411996*tau2 - 6.582724651790799797*delta3*tau2*tau - 10.17421240823160034*delta3*tau20*tau6 - 0.6652620261928799961*delta4 + 2.827822549952760056*delta4*tau2 + 4.930883534388300404*delta4*tau20*tau6 - 0.5671416341366160019*delta4*delta*tau2 - 0.8307233397384320428*delta4*delta2*tau20*tau6 + 9.794456308375400155*tau6 + 0.04170452661222047441*delta4*delta3*tau2 + 0.2326241138667191899*delta4*delta3*tau20*tau6 - 6.100523451393000229*tau6*tau + 0.007286832269659350436*delta4*delta4 - 0.004941628896799649291*delta4*delta4*delta*tau20*tau6 + 0.07884107375830799735*tau20*tau2 + 0.2511681684861600261*tau20*tau6)

def iapws97_dA_dtau_region3(tau, delta):
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta3*delta
    delta5 = delta4*delta
    
    tau11 = tau*tau # tau2
    tau3 = tau11*tau
    tau5 = tau3*tau11
    tau11 = tau5*tau5 # 10
    tau25 = tau11*tau11*tau5 # 25
    tau11 *= tau
    return (-2.530863095542799979*delta*tau - 6.914644684008599995*delta*tau5 + 13.27815659764770118*delta*tau11*tau3 - 10.91532008087319028*delta*tau11*tau5 - 15.37354157574320013*tau - 1.70429417648411996*tau*delta2 + 2.779959913892000145*tau*delta3 - 0.8895287085747799471*tau*delta4 + 0.1885215033301839982*tau*delta5*delta - 0.02700674448269600042*tau*delta5*delta2 + 0.001158459072561679893*tau*delta5*delta4 + 20.94439697430700065 + 29.38336892512619869*delta2*tau5 - 21.35183207987549991*delta2*tau5*tau + 0.867251811341387957*delta2*tau11*tau5*tau5 + 3.265186190320080506*delta2*tau25 - 8.075966009428000802*delta3*tau3 - 0.1314362194783408044*delta3*tau11*tau3*tau - 12.37496929107997978*delta3*tau25 + 3.622882828789319998*delta4*tau3 + 18.33583702287142003*delta4*tau25 + 0.1077051262633200029*delta5 - 0.9874086977686200584*delta5*tau*tau - 13.22647613070108008*delta5*tau25 + 4.273432396469860173*delta5*delta*tau25 + 18.33016345156779892*tau5*tau - 0.3856929791642719763*delta5*delta3*tau25 + 0.08400315222964860329*delta5*delta4*tau25 - 28.08078114861999808*tau5*tau3*tau - 0.0001655767979503699929*delta5*delta5 - 0.001168021375607189872*delta5*delta5*delta*tau25 + 14.4640436358203992*tau11 - 0.1945036694687546086*tau11*tau11)

def iapws97_d2A_d2tau_region3(tau, delta):
    delta2 = delta*delta
    delta3 = delta2*delta
    delta6 = delta3*delta3
    
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau10 = tau4*tau4*tau2
    tau24 = tau4*tau10*tau10
    return (-2.530863095542799979*delta - 34.57322342004299998*delta*tau4 + 185.8941923670678307*delta*tau10*tau2*tau - 174.6451212939710445*delta*tau10*tau4*tau - 1.974817395537240117*tau*delta3*delta2 - 15.37354157574320013 - 1.70429417648411996*delta2 + 146.9168446256310006*delta2*tau4 - 128.1109924792530137*delta2*tau4*tau + 18.21228803816914876*delta2*tau10*tau10 + 81.62965475800201887*delta2*tau24 + 2.779959913892000145*delta3 - 24.22789802828400241*delta3*tau2 - 1.97154329217511215*delta3*tau10*tau4 - 309.3742322769995212*delta3*tau24 - 0.8895287085747799471*delta3*delta + 10.86864848636795955*delta3*delta*tau2 + 458.395925571785483*delta3*delta*tau24 - 330.6619032675270091*delta3*delta2*tau24 + 109.9809807094067935*tau4*tau + 0.1885215033301839982*delta6 + 106.8358099117465088*delta6*tau24 - 0.02700674448269600042*delta6*delta - 9.642324479106799018*delta6*delta2*tau24 - 252.7270303375799756*tau4*tau4 + 0.001158459072561679893*delta6*delta3 + 2.100078805741214971*delta6*delta3*tau24 + 159.1044799940243877*tau10 - 0.02920053439017974636*delta6*delta3*delta2*tau24 - 4.279080728312601778*tau10*tau10*tau)


def iapws97_d2A_ddeltadtau_region3(tau, delta):
    tau2 = tau*tau
    delta2 = tau2*tau2 # abuse - tau4
    tau5 = delta2*tau
    tau14 = tau5*tau5 # 10
    tau25 = tau14*tau14*tau5
    tau14 *= delta2
    
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta3*delta
    delta5 = delta4*delta
    return (-3.408588352968239921*delta*tau + 58.76673785025239738*delta*tau5 - 42.70366415975099983*delta*tau5*tau + 1.734503622682775914*delta*tau14*tau5*tau2 + 6.530372380640161012*delta*tau25 - 2.530863095542799979*tau + 8.339879741676000435*tau*delta2 - 3.558114834299119789*tau*delta3 + 1.131129019981103934*tau*delta5 - 0.1890472113788720099*tau*delta5*delta + 0.0104261316530551186*tau*delta5*delta3 - 24.22789802828400241*delta2*tau2*tau - 0.3943086584350223855*delta2*tau14*tau - 37.12490787323994113*delta2*tau25 + 14.49153131515727999*delta3*tau2*tau + 73.34334809148568013*delta3*tau25 + 0.5385256313166000286*delta4 - 4.937043488843100292*delta4*tau2 - 66.13238065350540751*delta4*tau25 + 25.64059437881915926*delta5*tau25 - 6.914644684008599995*tau5 - 3.08554383331417581*delta5*delta2*tau25 + 0.7560283700668373186*delta5*delta3*tau25 - 0.001655767979503699984*delta5*delta4 - 0.01284823513167908729*delta5*delta5*tau25 + 13.27815659764770118*tau14 - 10.91532008087319028*tau14*tau2)
    
    
### Region 5
    
def iapws97_G0_region5(tau, pi):
    tau_inv = 1.0/tau
    return (tau*(6.85408416344340043 - 0.329616265389169993*tau)
            + tau_inv*(tau_inv*(0.369015349803330006 - 0.0248051489334660015*tau_inv) 
            - 3.1161318213925) + log(pi) - 13.1799836742010008)

def iapws97_dG0_dtau_region5(tau, pi):
    # does not depend on pi but leave as argument for consistency
    tau_inv = 1.0/tau
    return (-0.659232530778339987*tau + tau_inv*tau_inv*(tau_inv*(0.074415446800398008*tau_inv
            - 0.738030699606660012) + 3.1161318213925) + 6.85408416344340043)

def iapws97_d2G0_d2tau_region5(tau, pi):
    # does not depend on pi but leave as argument for consistency
    tau_inv = 1.0/tau
    return (tau_inv*tau_inv*tau_inv*(tau_inv*(2.21409209881998015 - 0.297661787201592032*tau_inv)
            - 6.232263642785) - 0.659232530778339987)

def iapws97_Gr_region5(tau, pi):
    tau3 = tau*tau*tau
    return (pi*(pi*tau3*(3.79194548229549995e-8*pi*tau*tau3 + (2.24400374094849992e-6 
            - 4.11632754534709986e-6*tau3*tau3)) + tau*(tau*(0.000901537616739440007 
            - 0.00502700776776479966*tau) + 0.00157364048552589993)))
    

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

def iapws97_d2Gr_d2pi_region5(tau, pi):
    tau2 = tau*tau
    return (tau2*tau*(tau2*tau2*(2.2751672893773001e-7*pi - 8.23265509069419973e-6*tau2)
                  + 4.48800748189699983e-6))

def iapws97_dGr_dtau_region5(tau, pi):
    tau2 = tau*tau
    return pi*(0.001803075233478880013*tau + 0.001573640485525899932
            + tau2*(-0.01508102330329439897 + pi*(6.732011222845499322e-6 
            + tau2*tau2*(-0.00003704694790812389877*tau2 
            + 2.654361837606849767e-7*pi))))

def iapws97_d2Gr_d2tau_region5(tau, pi):
    tau2 = tau*tau
# # 10 before
#    return (0.00180307523347888001*pi + tau*(0.0000134640224456909986*pi2 
#            - 0.0301620466065887979*pi + tau**4*pi2*(1.59261710256410975e-6*pi*
#                            - 0.00029637558326499119*tau2)))
    return pi*(0.001803075233478880013 + tau*(-0.03016204660658879794
            + pi*(0.00001346402244569099864 + tau2*tau2*(-0.0002963755832649911902*tau2
            + 1.592617102564109754e-6*pi))))

def iapws97_d2Gr_dpidtau_region5(tau, pi):
    tau2 = tau*tau
    return (0.001573640485525899932 + tau*(0.001803075233478880013
            + tau*(-0.01508102330329439897 + pi*(0.00001346402244569099864
            + tau2*tau2*(-0.00007409389581624779755*tau2 + 7.963085512820550889e-7*pi)))))


### Region 3 subregion density calls



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
        raise ValueError("For box (1,2,3,4) 273.15 K <= T <= 1073.15 K and P <= 100 MPa; "
                         "for box 5, 1073.15 K <= T <= 2273.15 K and P <= 50 MPa.")

def iapws97_region1_rho(T, P):
    # Useful to separate this out for confirming derivatives numerically
    pi = P*6.049606775559589e-08 #1/16.53E6
    tau = 1386.0/T
    dG_dpi = iapws97_dG_dpi_region1(tau, pi)
    return P/(R97*T*pi*dG_dpi)

def iapws97_region2_rho(T, P):
    pi = P*1e-6
    tau = 540.0/T
    dG_dpi = 1.0/pi + iapws97_dGr_dpi_region2(tau, pi) 
    return P/(R97*T*pi*dG_dpi)

def iapws97_region5_rho(T, P):
    pi = P*1e-6
    tau = 1000.0/T
    dG_dpi = 1.0/pi + iapws97_dGr_dpi_region5(tau, pi) 
    return P/(R97*T*pi*dG_dpi)


def iapws97_rho(T, P):
    r'''Calculate the density of water in kg/m^3 according to the IAPWS-97 
    standard.
    
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
    The range of validity of this formulation is as follows:
        
    For :math:`P \le 100 \text{ MPa}`:
        
    .. math::
        273.15 \text{ K} \le T \le 1073.15 \text{ K}

    For :math:`P \le 50 \text{ MPa}`:
        
    .. math::
        1073.15 \text{ K} \le T \le 2273.15 \text{ K}
        
    A ValueError is raised if the temperature or the pressure is out of bounds.

    IAPWS is implemented in four regions in the `T`-`P` domain:
    Region 1 (liquid), region 2 (gas and supercritical gas), region 5
    (high tempreature gas), and region 3 (near-critical).
    Significant discontinuities exist between the transitions of each regions.
    In region 3, there are 26 sub-regions and the correlation has the least 
    accuracy. 
    
    For many applications, the discontinuities in IF-97 can be problematic and
    the slower IAPWS-95 must be used. IAPWS-95 also has a wider range of
    applicability.
    
    Examples
    --------
    >>> iapws97_rho(648.6, 22.5e6)
    353.06081088726
    >>> iapws97_rho(330.0, 8e5)
    985.10498080770
    >>> iapws97_rho(823.0, 14e6)
    40.39293607288123
    >>> iapws97_rho(2000.0, 3e7)
    32.11456228328856
    
    References
    ----------
    .. [1] Cooper, JR, and RB Dooley. "Revised Release on the IAPWS Industrial 
       Formulation 1997 for the Thermodynamic Properties of Water and Steam." 
       The International Association for the Properties of Water and Steam 1 
       (2007): 48.
    '''
    region = iapws97_identify_region_TP(T, P)
    if region == 1:
        return iapws97_region1_rho(T, P)
    elif region == 2:
        return iapws97_region2_rho(T, P)
    elif region == 3:
        return iapws97_region3_rho(T, P)
    elif region == 5:
        return iapws97_region5_rho(T, P)
    else:
        raise ValueError("Out of bounds")


def iapws_97_Trho_err_region1(P, T, rho):
    pi_region1 = P*6.049606775559589e-08 #1/16.53E6
    tau_region1 = 1386.0/T
    dG_dpi_region1 = iapws97_dG_dpi_region1(tau_region1, pi_region1)
    rhol = P/(R97*T*pi_region1*dG_dpi_region1)
    err = rhol - rho
    d2G_dpi2_region1 = iapws97_d2G_dpi2_region1(tau_region1, pi_region1)
    derr = -d2G_dpi2_region1/(R97*T*dG_dpi_region1*dG_dpi_region1)
#    print(P, err, derr)
    return err, derr

def iapws_97_Trho_err_region2(P, T, rho):
    pi_region2 = P*1e-6
    tau_region2 = 540.0/T
    dG_dpi_region2 = 1/pi_region2 + iapws97_dGr_dpi_region2(tau_region2, pi_region2)
    rhog = P/(R97*T*pi_region2*dG_dpi_region2)
    err = rhog - rho

    d2G_dpi2_region2 = iapws97_d2Gr_d2pi_region2(tau_region2, pi_region2)
    d2G_dpi2_region2 -= 1e12/(P*P) # ideal part
    
    # checked numerically
    derr = -d2G_dpi2_region2/(R97*T*dG_dpi_region2*dG_dpi_region2)
#    print(P, T, rho, err, derr)
    return err, derr

def iapws_97_Trho_err_region5(P, T, rho):
    pi_region5 = P*1e-6
    tau_region5 = 1000.0/T
    dG_dpi_region5 = 1/pi_region5 + iapws97_dGr_dpi_region5(tau_region5, pi_region5)
    rhog = P/(R97*T*pi_region5*dG_dpi_region5)
    err = rhog - rho

    d2G_dpi2_region5 = iapws97_d2Gr_d2pi_region5(tau_region5, pi_region5)
    d2G_dpi2_region5 -= 1e12/(P*P) # ideal part
    derr = -d2G_dpi2_region5/(R97*T*dG_dpi_region5*dG_dpi_region5)
    return err, derr

def iapws97_P(T, rho):
    r'''Calculate the pressure of water according to the IAPWS-97 
    standard given a temperature `T` and mass density `rho`.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    rho : float
        Mass density of water, [kg/m^3]

    Returns
    -------
    P : float
        Pressure, [Pa]
        
    Notes
    -----    
    The range of validity of this formulation is as follows:
        
    For :math:`P \le 100 \text{ MPa}`:
        
    .. math::
        273.15 \text{ K} \le T \le 1073.15 \text{ K}

    For :math:`P \le 50 \text{ MPa}`:
        
    .. math::
        1073.15 \text{ K} \le T \le 2273.15 \text{ K}
        
    A ValueError is raised if the temperature or density is out of bounds.
    
    Newton's method with analytical derivatives is used here to solve these
    equations. The solver tolerance is as tight as it can be without causing
    wasted iterations that do not improve the result at all. Pressure changes
    quickly with density however, and some discrepancy between solvers is to be
    expected.
    
    For region 3, there are really two formulations present in IAPWS-97. There
    is a Helmholtz energy equation (Temperature and density dependent), and 
    also 26 separate backwards equations for `rho` which depent on `T` and `P`.
    The Helmholtz energy equation is much more accurate and does not have 
    discontinuities. The two sets of equations agree closely not not perfectly.
    By design, :obj:`iapws97_rho` implements the 26 T-P equations and this 
    implements the Helmholtz energy equation. This means that in region 3 
    solutions will not be consistent. For consistency requirements, IAPWS-95
    is recommended.
    
    This solver does not have any issues with multiple solutions. The solvers
    have been checked to achieve a relative solution tolerance of 5e-9 on 
    100 million points.
    
    Examples
    --------
    >>> iapws97_P(330.0, iapws97_rho(T=330.0, P=8e5))
    8e5
    >>> iapws97_P(823.0, 40.39293607288123)
    14e6
    >>> iapws97_P(T=2000.0, rho=32.11456228328856)
    3e7

    Region 3 point - does not implement the same equations as 
    :obj:`iapws97_rho`!
    
    >>> iapws97_P(648.6, iapws97_rho(T=648.6, P=22.5e6))
    22499974.093936257
    
    References
    ----------
    .. [1] Cooper, JR, and RB Dooley. "Revised Release on the IAPWS Industrial 
       Formulation 1997 for the Thermodynamic Properties of Water and Steam." 
       The International Association for the Properties of Water and Steam 1 
       (2007): 48.
    '''
    if T < 273.15:
        raise ValueError("T is under minimum value of 273.15 K")
    elif T <= 1073.15:
        if T <= 623.15:
            # region 2 to region 1 only - easy solver under 623K
            # Compute the density borders at the saturation pressure, and then decide which to pursue
            Psat = Psat_IAPWS(T)
            
            pi_region2 = Psat*1e-6
            tau_region2 = 540.0/T
            dG_dpi_region2 = 1.0/pi_region2 + iapws97_dGr_dpi_region2(tau_region2, pi_region2) 
            rhog_sat = Psat/(R97*T*pi_region2*dG_dpi_region2)

            pi_region1 = Psat*6.049606775559589e-08 #1/16.53E6
            tau_region1 = 1386.0/T
            dG_dpi_region1 = iapws97_dG_dpi_region1(tau_region1, pi_region1)
            rhol_sat = Psat/(R97*T*pi_region1*dG_dpi_region1)
            
            if rhog_sat < rho < rhol_sat:
                raise ValueError("Specified density is not a stable state at T")
            elif rho > rhol_sat:
                return newton(iapws_97_Trho_err_region1, Psat*10.0, fprime=True, bisection=True,
                              low=Psat, high=100e6, args=(T, rho), xtol=3e-12)                
            else:
                return newton(iapws_97_Trho_err_region2, Psat*.1, fprime=True, bisection=True,
                              low=Psat*1e-20, high=Psat, args=(T, rho), xtol=3e-12)
        P_region2_border = iapws97_boundary_2_3(T)
        pi_region2 = P_region2_border*1e-6
        tau_region2 = 540.0/T
        dG_dpi_region2 = 1.0/pi_region2 + iapws97_dGr_dpi_region2(tau_region2, pi_region2) 
        rhog_region2_border = P_region2_border/(R97*T*pi_region2*dG_dpi_region2)
        if rho < rhog_region2_border or P_region2_border > 100e6:
                return newton(iapws_97_Trho_err_region2, P_region2_border*.1, fprime=True, bisection=True,
                              low=P_region2_border*1e-20, high=P_region2_border, args=(T, rho), xtol=3e-12)
        else:
            # region 3
            tau = Tc/T
            delta = rho*rhoc_inv
            dA_ddelta = iapws97_dA_ddelta_region3(tau, delta)
            return dA_ddelta*delta*rho*R97*T
                
    elif T <= 2273.15:
        return newton(iapws_97_Trho_err_region5, 1e6, fprime=True, bisection=True,
                      low=1e-10, high=50e6, args=(T, rho), xtol=1e-12)
    else:
        raise ValueError("T is above maximum value of 2273.15 K")


def iapws_97_Prho_err_region1(T, P, rho):
    pi_region1 = P*6.049606775559589e-08 #1/16.53E6
    tau_region1 = 1386.0/T
    dG_dpi_region1 = iapws97_dG_dpi_region1(tau_region1, pi_region1)
        
    rhol = P/(R97*T*pi_region1*dG_dpi_region1)
    err = rhol - rho
    # what it is supposed to be
    drhol = (-16.53E6/(R97*T*T*dG_dpi_region1)
             + 16.53E6*1386.0*iapws97_d2G_dpidtau_region1(tau_region1, pi_region1)/(R97*T*T*T*dG_dpi_region1*dG_dpi_region1))
    return err, drhol

def iapws_97_Prho_err_region2(T, P, rho):
    pi_region2 = P*1e-6
    tau_region2 = 540.0/T
    dG_dpi_region2 = 1.0/pi_region2 + iapws97_dGr_dpi_region2(tau_region2, pi_region2)
        
    rhol = P/(R97*T*pi_region2*dG_dpi_region2)
    err = rhol - rho
    drhol = (-1e6/(R97*T*T*dG_dpi_region2)
             + 1E6*540.0*iapws97_d2Gr_dpidtau_region2(tau_region2, pi_region2)/(R97*T*T*T*dG_dpi_region2*dG_dpi_region2))
    return err, drhol

def iapws_97_Prho_err_region5(T, P, rho):
    pi_region5 = P*1e-6
    tau_region5 = 1000/T
    dG_dpi_region5 = 1.0/pi_region5 + iapws97_dGr_dpi_region5(tau_region5, pi_region5)
        
    rhol = P/(R97*T*pi_region5*dG_dpi_region5)
    err = rhol - rho
    drhol = (-1e6/(R97*T*T*dG_dpi_region5)
             + 1E6*1000*iapws97_d2Gr_dpidtau_region5(tau_region5, pi_region5)/(R97*T*T*T*dG_dpi_region5*dG_dpi_region5))
    return err, drhol


def iapws_97_Prho_err_region3(T, P, rho):
    tau = Tc/T
    delta = rho*rhoc_inv
    dA_ddelta = iapws97_dA_ddelta_region3(tau, delta)
    P_calc = dA_ddelta*delta*rho*R97*T
    err = P_calc - P

    d2A_ddeltadtau = iapws97_d2A_ddeltadtau_region3(tau, delta)
    
    derr = R97*rho**2*dA_ddelta/rhoc - R97*Tc*rho**2*d2A_ddeltadtau/(T*rhoc)
    return err, derr

def iapws97_T(P, rho):
    r'''Calculate the temperature of water according to the IAPWS-97 
    standard given a pressure `P` and mass density `rho`.
    
    Parameters
    ----------
    P : float
        Pressure, [Pa]
    rho : float
        Mass density of water, [kg/m^3]

    Returns
    -------
    T : float
        Temperature, [K]
        
    Notes
    -----    
    The range of validity of this formulation is as follows:
        
    For :math:`P \le 100 \text{ MPa}`:
        
    .. math::
        273.15 \text{ K} \le T \le 1073.15 \text{ K}

    For :math:`P \le 50 \text{ MPa}`:
        
    .. math::
        1073.15 \text{ K} \le T \le 2273.15 \text{ K}
        
    A ValueError is raised if the pressure or density is out of bounds.
    
    Newton's method with analytical derivatives is used here to solve these
    equations. The solver tolerance is as tight as it can be without causing
    wasted iterations that do not improve the result at all. 
    
    Due to water's unique density curve, there is a temperature region
    spaning 273.15 K to 280.005 K where there are two solutions. No guarantee
    is made as to which solution will be returned.
    
    Examples
    --------
    >>> iapws97_T(8e5, iapws97_rho(T=330.0, P=8e5))
    330.0
    >>> iapws97_T(14e6, 40.39293607288123)
    823.0
    >>> iapws97_T(P=3e7, rho=32.11456228328856)
    2000.0
    
    References
    ----------
    .. [1] Cooper, JR, and RB Dooley. "Revised Release on the IAPWS Industrial 
       Formulation 1997 for the Thermodynamic Properties of Water and Steam." 
       The International Association for the Properties of Water and Steam 1 
       (2007): 48.
    '''
    solve_region = 0
    if P > 100e6:
        raise ValueError("P is above maximum value of 100 MPa")

    if P < 50E6:
        # Calculate the 2-5 border using region 5's equations
        pi_5_border = P*1e-6
        tau_5_border = 1000.0/1073.15
        dG_dpi_5_border = 1.0/pi_5_border + iapws97_dGr_dpi_region5(tau_5_border, pi_5_border) 
        rho_25_border = P/(R97*1073.15*pi_5_border*dG_dpi_5_border)

        # get rho minimum at 1073.15
        if rho <= rho_25_border:
            rho_5end_border = iapws97_rho(2273.15, P)
            if rho < rho_5end_border:
                raise ValueError("Density is lower then region 5 limit")
            else:
                solve_region = 5
    if solve_region == 0:
        no_region_3 = False
        if P > 13918839.778885445:
            T_23 = iapws97_boundary_2_3_reverse(P)
            if T_23 < 623.15:
                no_region_3 = True
        else:
            no_region_3 = True
        if no_region_3:
            if P < 0.005706860511498897: # where the vapor pressure equation dies
                solve_region = 2 # always gas
                T_23 = 273.15
            else:
                Tsat = Tsat_IAPWS(P)
                if Tsat < 273.15:
                    solve_region = 2  # always gas
                    T_23 = 273.15
                else:
                    pi_region1_sat = P*6.049606775559589e-08 #1/16.53E6
                    tau_region1_sat = 1386.0/Tsat
                    dG_dpi_region1_sat = iapws97_dG_dpi_region1(tau_region1_sat, pi_region1_sat)
                    rho1_sat = P/(R97*Tsat*pi_region1_sat*dG_dpi_region1_sat)

                    pi_region2_sat = P*1e-6
                    tau_region2_sat = 540.0/Tsat
                    dG_dpi_region2_sat = 1.0/pi_region2_sat + iapws97_dGr_dpi_region2(tau_region2_sat, pi_region2_sat) 
                    rho2_sat = P/(R97*Tsat*pi_region2_sat*dG_dpi_region2_sat)

#                     if rho2_sat < rho < rho1_sat:
#                         raise ValueError("At specified pressure, density is not a stable solution")
                    if rho > rho2_sat:
                        solve_region = 1
                    else:
                        solve_region = 2
                        T_23 = Tsat # for initial guess bondary only

        else:
            rho_23_side3 = iapws97_rho(T_23 * (1 + 1e-12), P)
            
            # Calculate the 2-5 border using region 2's equations
            pi_region2_25_on2 = P*1e-6
            tau_region2_25_on2 = 540.0/1073.15
            dG_dpi_region2_25_on2 = 1.0/pi_region2_25_on2 + iapws97_dGr_dpi_region2(tau_region2_25_on2, pi_region2_25_on2) 
            rho2_25_on2 = P/(R97*1073.15*pi_region2_25_on2*dG_dpi_region2_25_on2)
        
            if rho2_25_on2 <= rho <= rho_23_side3:
                solve_region = 2
    if solve_region == 0:
        rho_13 = iapws97_rho(623.15-1e-12, P)
        if rho_23_side3 <= rho <= rho_13:
            solve_region = 3
        if solve_region == 0:
            rho_Tmin_border = iapws97_rho(273.15, P)
            if rho_13 <= rho <= rho_Tmin_border:
                solve_region = 1
            else:
                # Sometimes liquid water has a higher density at a higher temperature. Weird!
                solve_region = 1
    if solve_region == 5:
        return newton(iapws_97_Prho_err_region5, 1673.15, fprime=True, bisection=True,
                      low=1073.15, high=2273.15, args=(P, rho), xtol=1e-12)
    elif solve_region == 2:
        return newton(iapws_97_Prho_err_region2, 0.5*(T_23 + 1073.15), fprime=True, bisection=True,
                      low=T_23 , high=1073.15, args=(P, rho), xtol=1e-12)
    elif solve_region == 3:
        return newton(iapws_97_Prho_err_region3, 0.5*(T_23 + 623.15), fprime=True, bisection=True,
                      low=623.15 , high=T_23, args=(P, rho), xtol=1e-12)
    elif solve_region == 1:
        try:
            Tsat = Tsat_IAPWS(P)
        except:
            Tsat = 623.15
        return newton(iapws_97_Prho_err_region1, 0.5*(Tsat + 273.15), fprime=True, bisection=True,
                      low=273.15 , high=Tsat, args=(P, rho), xtol=1e-12)
    else:
        raise ValueError("Could not detect region")

### IAPWS95

### IAPWS 95 Initial Guesses

def rhol_sat_IAPWS(T):
    r'''Calculates saturation liquid mass density of water using the IAPWS
    SR1-86(1992) [1]_ [2]_ explicit equation.

    .. math::
        \frac{\rho^{sat}_l}{\rho_c} = 1 + b_1\tau^{1/3} + b_2\tau^{2/3}
        + b_3 \tau^{5/3} + b_4\tau^{16/3} + b_5\tau^{43/3} + b_6\tau^{110/3}
        
    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature of water, [K]

    Returns
    -------
    rhol_sat : float
        Saturation liquid mass density of water [kg/m^3]

    Notes
    -----
    This equation is fit to experimental data to within its accuracy. It does 
    not satisfy the equilibrium conditions for the IAPWS-95 or IAPWS-97 
    formulations.
    
    The values of the constants are as follows: 
        
    b1 = 1.99274064; 
    b2 = 1.09965342; 
    b3 = -0.510839303; 
    b4 = -1.75493479; 
    b5 = -45.5170352; 
    b6 = -6.74694450e5


    Examples
    --------
    >>> rhol_sat_IAPWS(300.)
    996.5089712803

    References
    ----------
    .. [1] IAPWS, Secretariat, B Dooley, and EPRI. "Revised Supplementary 
       Release on Saturation Properties of Ordinary Water Substance", 1992.
    .. [2] Wagner, Wolfgang, and A. Pruss. "International Equations for the
       Saturation Properties of Ordinary Water Substance. Revised According to 
       the International Temperature Scale of 1990. Addendum to J. Phys. Chem.
       Ref. Data 16, 893 (1987)." Journal of Physical and Chemical Reference 
       Data 22, no. 3 (May 1, 1993): 783-87. https://doi.org/10.1063/1.555926.
    '''
    tau = 1.0 - T*Tc_inv
    
    tau_cbrt = tau**(1.0/3.0)
    
    ratio = 1.0 + 1.99274064*tau_cbrt
    tau_cbrt4 = tau_cbrt*tau_cbrt # still 2 for first term
    
    ratio += 1.09965342*tau_cbrt4 # still b2*tau^(2/3)
    
    tau_cbrt4 = tau_cbrt4*tau_cbrt4
    
    ratio += -0.510839303*tau_cbrt4*tau_cbrt
    
    tau_cbrt8 = tau_cbrt4*tau_cbrt4
    
    ratio += -1.75493479*tau_cbrt8*tau_cbrt8
    
    tau_cbrt4 = tau_cbrt8*tau_cbrt # repurpse 4 as 9
    
    tau_cbrt = tau_cbrt4*tau_cbrt8 # 17 - repurpose variable as tau_cbrt34
    tau_cbrt *= tau_cbrt # 17+17 = 34
    
    ratio += tau_cbrt*(-45.5170352*tau_cbrt4 -6.74694450e5*tau_cbrt*tau_cbrt*tau_cbrt8)
#    ratio += -6.74694450e5*tau_cbrt*tau_cbrt*tau_cbrt*tau_cbrt8
    return ratio*rhoc

def rhog_sat_IAPWS(T):
    r'''Calculates saturation vapor mass density of water using the IAPWS
    SR1-86(1992) [1]_ [2]_ explicit equation.

    .. math::
        \log \left(\frac{\rho^{sat}_g}{\rho_c}\right) = 1 + c_1\tau^{2/6} + c_2\tau^{4/6}
        + c_3 \tau^{8/6} + c_4\tau^{18/6} + c_5\tau^{37/6} + c_6\tau^{71/6}
        
    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature of water, [K]

    Returns
    -------
    rhog_sat : float
        Saturation vapor mass density of water [kg/m^3]

    Notes
    -----
    This equation is fit to experimental data to within its accuracy. It does 
    not satisfy the equilibrium conditions for the IAPWS-95 or IAPWS-97 
    formulations.
    
    The values of the constants are as follows: 
        
    c1 = -2.03150240;
    c2 = -2.68302940;
    c3 = -5.38626492;
    c4 = -17.2991605;
    c5 = -44.7586581;
    c6 = -63.9201063


    Examples
    --------
    >>> rhog_sat_IAPWS(300.)
    0.0255887212886

    References
    ----------
    .. [1] IAPWS, Secretariat, B Dooley, and EPRI. "Revised Supplementary 
       Release on Saturation Properties of Ordinary Water Substance", 1992.
    .. [2] Wagner, Wolfgang, and A. Pruss. "International Equations for the
       Saturation Properties of Ordinary Water Substance. Revised According to 
       the International Temperature Scale of 1990. Addendum to J. Phys. Chem.
       Ref. Data 16, 893 (1987)." Journal of Physical and Chemical Reference 
       Data 22, no. 3 (May 1, 1993): 783-87. https://doi.org/10.1063/1.555926.
    '''
    tau = 1.0 - T*Tc_inv
    
    tau_6rt = tau**(1.0/6.0)
    tau_6rt2 = tau_6rt*tau_6rt
    ratio = -2.03150240*tau_6rt2
    
    tau_6rt8 = tau_6rt2*tau_6rt2 # start it off as 4
    ratio += -2.68302940*tau_6rt8
    
    
    tau_6rt8 *= tau_6rt8
    ratio += -5.38626492*tau_6rt8
    
    tau_6rt16 = tau_6rt8*tau_6rt8
    tau_6rt18 = tau_6rt16*tau_6rt2 
    
    ratio += -17.2991605*tau_6rt18
    
    tau_6rt2 = tau_6rt18*tau_6rt18*tau_6rt # 37 - reuse tau_6rt2
    ratio += tau_6rt2*(-44.7586581 - 63.9201063*tau_6rt16*tau_6rt18) # 71
    
    return exp(ratio)*rhoc

### IAPWS 95 fundamental derivatives

def iapws95_A0(tau, delta):
    r'''Calculates the ideal gas Helmholtz energy of water according to the
    IAPWS-95 standard.
    
    .. math::
        \phi^\circ = \ln \delta + n_1 + n_2\tau + n_3\ln \tau
        + \sum_{i=4}^8 n_i \ln \left[1 - \exp(-\gamma_i \tau)   \right]

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    A0 : float
        Ideal gas dimensionless Helmholtz energy A/(RT) [-]

    Notes
    -----
    This implementation is checked to have a relative error always under 1e-15.
            
    Examples
    --------
    >>> iapws95_A0(647.096/300.0, 999.0/322)
    9.537075529761053
    '''
    # only way to optimize is likely to save transcendentals between calls
    # It should also be possible to replace the tau bit with polynomial fits 
    # which are nearly bit for bit identical, including the derivatives.
    return (6.68321052759320011*tau + log(delta) + 3.00632*log(tau) 
            + 0.24873*log(1. - exp(-27.5075105*tau)) 
            + 0.96956*log(1. - exp(-9.24437796*tau)) 
            + 1.2795*log(1. - exp(-7.74073708*tau)) 
            + 0.97315*log(1. - exp(-3.53734222*tau))
            + 0.012436*log(1. - exp(-1.28728967*tau))
            - 8.32044648374970031)

def iapws95_dA0_dtau(tau, delta):
    r'''Calculates the first derivative of ideal gas Helmholtz energy of water
    with respect to `delta` according to the IAPWS-95 standard.
    
    .. math::
        \frac{\partial \phi^\circ}{\partial \tau} = n_2 + \frac{n_3}{\tau}
        + \sum_{i=4}^8 n_i\gamma_i \left[\left(1-\exp(-\gamma_i \tau)
        \right)^{-1} - 1 \right]

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    dA0_dtau : float
        First derivative of ideal gas dimensionless Helmholtz energy A/(RT)
        with respect to `tau` [-]

    Notes
    -----
    This implementation is checked to have a relative error always under 1e-15.
            
    Examples
    --------
    >>> iapws95_dA0_dtau(647.096/300.0, 999.0/322)
    8.079705548882
    '''
    return (-22.4843580635585205 
            + 0.01600873433612/(1.0 - exp(-1.28728967*tau))
            + 3.44236458139299994/(1.0 - exp(-3.53734222*tau)) 
            + 9.90427309386000054/(1.0 - exp(-7.74073708*tau)) 
            + 8.96297909489759981/(1.0 - exp(-9.24437796*tau)) 
            + 6.84194308666500017/(1.0 - exp(-27.5075105*tau))
            + 3.00632/tau)

def iapws95_d2A0_dtau2(tau, delta):
    r'''Calculates the second derivative of ideal gas Helmholtz energy of water
    with respect to `delta` according to the IAPWS-95 standard.
    
    .. math::
        \frac{\partial^2 \phi^\circ}{\partial \tau^2} = \frac{n_3}{\tau^2}
        + \sum_{i=4}^8 n_i\gamma_i ^2 \exp(-\gamma_i \tau) 
        \left[\left(1-\exp(-\gamma_i \tau)
        \right)^{-2}\right]

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    d2A0_dtau2 : float
        Second derivative of ideal gas dimensionless Helmholtz energy A/(RT)
        with respect to `tau` [-]

    Notes
    -----
    This implementation is checked to have a relative error always under 1e-15.
            
    Examples
    --------
    >>> iapws95_d2A0_dtau2(647.096/300.0, 999.0/322)
    -0.653543047751809
    '''
    x0 = exp(-27.5075105*tau)
    x1 = exp(-9.24437796*tau)
    x2 = exp(-7.74073708*tau)
    x3 = exp(-3.53734222*tau)
    x4 = exp(-1.28728967*tau)
    return (-188.204821296839867*x0/((1.0 - x0)*(1.0 - x0))
            - 82.8571664008121047*x1/((1.0 - x1)*(1.0 - x1))
            - 76.6663739880884236*x2/((1.0 - x2)*(1.0 - x2))
            - 12.1768215703940861*x3/((1.0 - x3)*(1.0 - x3))
            - 0.0206078783406615819*x4/((1.0 - x4)*(1.0 - x4))
            - 3.00632/(tau*tau))

def iapws95_d3A0_dtau3(tau, delta):
    r'''Calculates the third derivative of ideal gas Helmholtz energy of water
    with respect to `delta` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    d3A0_dtau3 : float
        Third derivative of ideal gas dimensionless Helmholtz energy A/(RT)
        with respect to `tau` [-]

    Notes
    -----
    This implementation is checked to have a relative error always under 1e-15.
    This equation is not explicitly in IAPWS-95, but is needed to compute some
    second derivatives.
            
    Examples
    --------
    >>> iapws95_d3A0_dtau3(647.096/300.0, 999.0/322)
    0.6222542507278
    '''
    tot = 0.0
    x0 = exp(-27.5075105*tau)
    x1 = 1.0/(1.0 - x0)
    x2 = exp(-9.24437796*tau)
    x3 = 1.0/(1.0 - x2)
    x4 = exp(-7.74073708*tau)
    x5 = 1.0/(1.0 - x4)
    x6 = exp(-3.53734222*tau)
    x7 = 1.0/(1.0 - x6)
    x8 = exp(-1.28728967*tau)
    x9 = 1.0/(1.0 - x8)
    return (x1*x1*x0*(5177.04609797344619
            + 10354.0921959468924*x0*x1)
            + (x3*x3*x2)*(765.962962903719927
            + 1531.92592580743985*x2*(x3))
            + x5*x5*x4*(593.454243918743487
            + 1186.90848783748697*x4*(x5))
            + x7*x7*x6*(86.1471700927234139*x6*(x7)
            + 43.0735850463617069)
            + x9*x9*x8*(0.0265283089085503951 
            + 0.0530566178171007902*x8*(x9))
            + 6.01264/(tau*tau*tau))



def iapws95_Ar(tau, delta):
    r'''Calculates the residual Helmholtz energy of water according to the
    IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    Ar : float
        Residual Helmholtz energy A/(RT) [-]

    Notes
    -----
    This is an optimized implementatation taking 9 exp calls, 4 sqrts, 
    and 3 powers. It was generated using SymPy's CSE functionality, with 
    select polynomial optimizations by hand as well. It is over
    10x faster than a naive implementation.
    
    This implementation has been tested against a straightforward 
    implementation with the equations given in IAPWS-95.
    
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-10
    kg/m^3 to 5000 kg/m^3, 4E6 points were evaluated. The mean relative
    error was 5.0416E-15, with a maximum relative error of 1.118E-9 and a 
    standard deviation of 5.773e-13.    

    Over the same range, the model was evaluated to a precision of 50 
    decimal places with `mpmath`, and on 90000 points, the mean relative
    error was 3.14E-15, with a maximum relative error of 3.54e-12 and a 
    standard deviation of 3.017E-14.

    This comparison indicates that this implementation is more accurate than 
    the straightforward implementation.
            
    Examples
    --------
    >>> iapws95_Ar(647.096/300.0, 999.0/322)
    -9.57577716026768
    '''
    _sqrt, _exp = sqrt, exp
    taurt = _sqrt(tau)
    tau_quarter = _sqrt(taurt)
    tau_eighth = _sqrt(tau_quarter) # tau checked, is not causing the small discrepancies.
    deltam1sqr = (delta - 1.0)
    deltam1sqr *= deltam1sqr
    taum1sqr = (tau - 1.0)
    taum1sqr *= taum1sqr
    tau2 = tau*tau
    tau3 = tau*tau2
    tau4 = tau*tau3
    tau6 = tau4*tau2
    x29 = tau4*tau4
    tau9 = tau6*tau3
    tau10 = tau*tau9
    tau11 = tau10*tau
    tau13 = tau10*tau3
    tau50 = tau10*tau10
    tau50 *= tau50*tau10
    delta2 = delta*delta
    delta4 = delta2*delta2
    delta5 = delta*delta4
    delta6 = delta*delta5
    delta8 = delta6*delta2
    delta9 = delta*delta8
    delta3 = delta*delta2
    expnd = _exp(-delta)
    expnd2 = _exp(-delta2)
    expnd3 = _exp(-delta3)
    expnd6 = _exp(-delta6)
    tauexpnd = tau*expnd
    tau4expnd = expnd*tau4
    deltaexpnd2 = delta*expnd2
    tau7expnd2 = tau6*tau*expnd2
    delta3expnd = delta3*expnd
    delta8expnd2 = delta8*expnd2
    tau10expnd2 = tau10*expnd2
    tau23expnd3 = tau10*tau13*expnd3
    delta5expnd6 = expnd6*delta5
    delta9expnd2 = delta9*expnd2
    x32 = -20.0*deltam1sqr
    x33 = (0.826446280991736*tau - 1.0)
    x33 = delta2*_exp(x32 - 219.615*x33*x33)
    x50 = (-tau + 0.32*(deltam1sqr)**1.66666666666666674 + 1.0)
    x51 = 0.8*tau - 1.0
    x35 = 0.2*_sqrt(deltam1sqr)*deltam1sqr*deltam1sqr*deltam1sqr + x50*x50
    x35_n05 = x35**-0.05
    return (delta*(delta*(delta9*(delta*(delta*(delta*(-6.26395869124539993e-10*delta*tauexpnd 
                                - 0.000062689710414685001*tau10*_exp(-delta4))
           - 1.32511800746680002e-12*tau13*expnd) - 0.0000163885683425300002*x29*expnd2)
           + 3.65821651442040008e-7*tau4expnd)
           
           + taurt*(-0.26145533859358*tau_quarter + 0.318025093454180008)
           -tauexpnd*(0.25709043003437998*tau4 + 0.192327211560020001)
           - 0.00781997516879810034*tau_quarter*tau_eighth*delta
           )
           + 7.89576347228280007*tau**0.875
           
           + tau11*(tau11*tau11*tau11*delta5expnd6*(0.317774973307380026*tau2
           - 0.199057183544080002)
           + expnd*(1.15379964229510002e-9*delta9
           - 0.0000662126050396869941*tau))
           
           # These two terms catastrophic cancel somewhat; cannot change them any at all
           + 0.0349940054637650003*tau11*tau11*delta3*expnd3
           + 0.0436136157238109987*tau4*tau*tau11*delta2*expnd3
           
           + tau50*(expnd6*delta2*(-0.118411824259810006*delta3
            - 5.57111185656449973e-10))
           + tau*(0.00880894931021340005*delta3
           + 0.00833265048807130086*delta8expnd2
           + 0.0176114910087519991*deltaexpnd2
           + 31.5461402377809996*x33
           - 8.78032033035609949)
           
           +delta3*(expnd2*tau3*(tau4*(-0.0313587007125490022
                                 - 0.743159297103409999*tau3)
               + 0.0049969146990805997)
           # tau23expnd3 needs to do in both
           - 0.0767881978446210006*tau23expnd3
            +delta*(0.0224462773320059997*tau23expnd3
           - 7.59413770881439999e-6*tau9*expnd
           + 0.478073299154800013*tau10expnd2))
           
            - 0.107936009089319995*tau7expnd2
           + 0.000158703083241569997*tau9*delta9expnd2 
           + 0.221322951675460011*tau9*deltaexpnd2
           - 0.40247669763527999*tau10*deltaexpnd2 
           - 0.0400928289258069975*tau2*delta3expnd 
           - 0.029052336009585001*tau2*delta8expnd2
           + 3.93434226032540015e-7*delta3expnd*tau13
           + 0.000562509793518880044*delta6*tau3*expnd
           + 0.0141806344006170006*delta6*tau10expnd2
           + 0.0386150855742060026*tau3*delta8expnd2
           - 0.0000156086522571349985*delta8*tau4expnd 
           + 0.580833999857589989*delta2*tau10expnd2
           - 2521.31543416949989*delta2*tau4*_exp(x32 - 390.625*x51*x51)
           + 0.160748684862510011*delta2*tau4expnd
           - 0.0203934865137039983*delta8expnd2*tau4
           - 0.136364351103430009*tau10expnd2*delta5
            + 0.0205279408959480013*delta5*tau6*expnd2 
            + 0.204338109509650007*expnd*tau6
           + 0.00199555719795409996*delta9expnd2*tau6
           
           - 0.00165540500637340006*delta8expnd2*x29 
           - 31.3062603234350014*x33
           - 0.148746408567240002*x35*x35_n05*x35_n05*x35_n05*_exp(-28.0*deltam1sqr - 700.0*taum1sqr)
           + 0.318061108784439994*x35*x35_n05*_exp(-32.0*deltam1sqr - 800.0*taum1sqr)
           - 0.668565723079650009*tau4expnd + 0.0125335479355230001/taurt))


def iapws95_dAr_ddelta(tau, delta):
    r'''Calculates the first derivative of residual Helmholtz energy of water
    with respect to `delta` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    dAr_ddelta : float
        First derivative of residual Helmholtz energy A/(RT) with respect to
        `delta`, [-]

    Notes
    -----
    This is an optimized implementatation taking 11 exp calls, 4 sqrts, 
    and 3 powers. It was generated using SymPy's CSE functionality, with 
    select polynomial optimizations by hand as well. It is over
    10x faster than a naive implementation.
    
    This implementation has been tested against a straightforward 
    implementation with the equations given in IAPWS-95.
    
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-10
    kg/m^3 to 5000 kg/m^3, 4E6 points were evaluated. The mean relative
    error was 4.033E-15, with a maximum relative error of 3.8765e-10 and a 
    standard deviation of 3.189e-13.    
    
    Over the same range, the model was evaluated to a precision of 50 
    decimal places with `mpmath`, and on 90000 points, the mean relative
    error was 6.046E-15, with a maximum relative error of 3.39E-10 and a 
    standard deviation of 7.056E-13.
    
    Examples
    --------
    >>> iapws95_dAr_ddelta(647.096/300.0, 999.0/322)
    -0.3093321202374
    '''
    _sqrt, _exp = sqrt, exp
    taurt = _sqrt(tau)
    tau_quarter = _sqrt(taurt)
    tau_eighth = _sqrt(tau_quarter) # tau checked, is not causing the small discrepancies.
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta*delta3
    delta5 = delta*delta4
    delta6 = delta*delta5
    delta8 = delta6*delta2
    delta9 = delta*delta8
#    10, 11, 12, 13, 14 also used for delta but not needed
    x6 = delta - 1.0
    x7 = _exp(-delta)
    x25 = _exp(-delta2)
    x37 = _exp(-delta6)
    x32 = _exp(-delta3)
    x8 = x6*x7
    tau2 = tau*tau
    tau3 = tau*tau2
    tau4 = tau*tau3
    tau6 = tau3*tau3
    tau7 = tau6*tau
    tau9 = tau6*tau3
    tau10 = tau*tau9
    tau12 = tau10*tau2
    tau13 = tau10*tau3
    tau23 = tau10*tau13
    x10 = tau*x7
    x11 = delta*(delta - 2.0)
    x12 = tau4*x7
    x14 = delta3*x7*(delta - 4.0)
    x21 = delta8*tau4
    x24 = delta2 + delta2
    x26 = tau*tau6*x25
    x27 = delta*(delta2 - 1.0)
    x28 = x25*x27
    x30 = x25*tau10
    x33 = delta3*(delta2 - 2.0)
    x35 = delta5*(delta2 - 3.0)
    x38 = delta5*(delta6 - 1.0)
    x39 = x37*x38
    x40 = tau23*tau23*tau4*x37
    x41 = delta9*x25*(delta2 - 5.0)
    x42 = tau4*tau4*x25
    x43 = 3.0*delta3
    x44 = delta3*(x43 - 4.0)
    x45 = tau23*x32
    x46 = x24 - 9.0
    x47 = x25*x46
    x48 = delta8*x47
    x50 = delta - 1.0
    x51 = x50*x50
    dn1_2_23 = x51**(2.0/3.0)
    x65 = x51*dn1_2_23
    x52 = -20.0*x51
    x53 = delta3*(-40.0*delta + 40.0 + 3.0/delta)
    c54 = (tau - 1.21)
    x54 = x53*_exp(x52 - 150.0*c54*c54)
    x51_2_x51sqrt = x51*x51*_sqrt(x51)
    y50 = 0.32*x65 + 1.0 - tau
    x56 = 0.2*x51*x51_2_x51sqrt + y50*y50
    x57 = 32.0*x51
    x58 = (tau - 1.0)
    x58 *= x58
    x59 = 800.0*x58
    x62 = delta*x6
    x63 = 28.0*x51
    x64 = 700.0*x58
    c51 = _exp(-x57 - x59)
    c52 = _exp(-x63 - x64)
    x100 = (tau - 1.25)
    x101 = x56**(-0.05)
    x67 = x62*(1.4*x51_2_x51sqrt + (-2.13333333333333333*tau + 0.68266666666666666*x65 + 2.13333333333333333)*dn1_2_23)
    return (-3.6582165144204e-7*delta9*delta*x12*(delta - 11.0) + 3.277713668506e-5*delta9*delta2*x42*(delta2 - 6)
            + 1.3251180074668e-12*delta9*delta3*tau13*x7*(delta - 13.0)
            + 0.00012537942082937*delta9*delta4*tau10*(delta4 + delta4 - 7.0)*_exp(-delta4)
            + 6.2639586912454e-10*delta9*delta5*x10*(delta - 15.0) - 0.023459925506394301*tau_quarter*tau_eighth*delta2
            - 0.52291067718716*tau_quarter*taurt*delta 
            + 7.89576347228280007*tau**0.875#tau_quarter*taurt*tau_eighth # improves average but increases max.
            + 0.25709043003438*tau4*tau*x11*x7
            - 1.1537996422951e-9*tau*tau10*delta9*x7*(delta - 10.0) + 6.6212605039687e-5*tau12*x8
            - 0.130840847171432989*tau9*tau7*delta2*x32*(delta3 - 1.0) - 0.0349940054637650003*tau10*tau12*x32*x44
            + 1.19434310126448007*tau23*tau12*tau9*x39 - 1.90664983984428016*tau23*tau23*x39 + 0.636050186908360016*taurt*delta
            - 0.0352229820175039982*tau*x28 + 0.0352357972408536002*tau*delta3 - 0.00833265048807130086*tau*x48
            + 31.546140237781*tau*x54 - 8.7803203303561*tau + 0.19232721156002*x10*x11
            - 0.16074868486251*x12*delta2*(delta - 3.0) + 0.0400928289258069975*tau2*x14 + 0.029052336009585001*tau2*x48
            - 3.9343422603254e-7*x14*tau13 + 7.5941377088144e-6*delta4*tau9*x7*(delta - 5.0)
            - 0.4780732991548*delta4*x30*(x24 - 5.0) - 0.0224462773320059997*delta4*x45*(x43 - 5.0)
            - 0.442645903350920022*tau9*x28 - 0.000317406166483139994*tau9*x41 - 0.000562509793518880044*delta6*tau3*x7*(delta - 7.0)
            - 0.0141806344006170006*delta6*x30*(x24 - 7.0) - 0.0099938293981611994*tau3*x25*x33
            - 0.0386150855742060026*tau3*x48 + 0.00165540500637340006*delta8*x42*x46 + 0.0203934865137039983*x21*x47
            + 1.5608652257135e-5*x21*x7*(delta - 9.0) - 0.0410558817918960026*x25*x35*tau6 + 0.0627174014250980044*x26*x33
            + 0.10793600908932*x26*(x24 - 1.0) + 0.804953395270559979*x27*x30 - 0.58083399985759*delta2*x30*(x24 - 3.0)
            + 1.67133355696935002e-9*delta2*x40*(delta6 + delta6 - 1.0) + 1.48631859420682*x30*x33 + 0.272728702206860019*x30*x35
            + 0.710470945558860034*x38*x40 - 0.00399111439590819992*x41*tau6 + 0.0767881978446210006*x44*x45
            - 2521.3154341695*tau4*x53*_exp(x52 - 250.0*x100*x100) + 0.668565723079650009*tau4*x8
            - 31.3062603234350014*x54 - 0.14874640856724*x56*x101*x101*x101*(-56.0*x62*_exp(-28.0*x51 - 700.0*x58)
            + c52) + 0.31806110878444*x56*x101*(-64.0*x62*_exp(-32.0*x51 - 800.0*x58)
            + c51) - 0.126434447282154*x101*x101*x101*x67*c52
            + 0.302158053345218003*x101*x67*c51 - 0.204338109509650007*x8*tau6 + 0.0125335479355230001/taurt)

            
            

def iapws95_d2Ar_ddelta2(tau, delta):
    r'''Calculates the second derivative of residual Helmholtz energy of water
    with respect to `delta` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    d2Ar_ddelta2 : float
        Second derivative of residual Helmholtz energy A/(RT) with respect to
        `delta`, [-]

    Notes
    -----
    This is an optimized implementatation taking 4 exp calls, 4 sqrts, 
    and 2 powers. It was generated using SymPy's CSE functionality, with 
    select polynomial optimizations by hand as well. It is over
    10x faster than a naive implementation.
    
    This implementation has been tested against a straightforward 
    implementation with the equations given in IAPWS-95.
    
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-10
    kg/m^3 to 5000 kg/m^3, 4E6 points were evaluated. The mean relative
    error was 9.566e-16, with a maximum relative error of 1.0518E-10 and a 
    standard deviation of 6.20265E-14.    
    
    Over the same range, the model was evaluated to a precision of 50 
    decimal places with `mpmath`, and on 250000 points, the mean relative
    error was 1.039E-15, with a maximum relative error of 2.431E-11 and a 
    standard deviation of 5.31708E-14.

    Examples
    --------
    >>> iapws95_d2Ar_ddelta2(647.096/300.0, 999.0/322)
    1.7862535141735987
    '''
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

def iapws95_d3Ar_ddelta3(tau, delta):
    r'''Calculates the third derivative of residual Helmholtz energy of water
    with respect to `delta` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    d3Ar_ddelta3 : float
        Third derivative of residual Helmholtz energy A/(RT) with respect to
        `delta`, [-]

    Notes
    -----
    This is an optimized implementatation.
    It was generated using SymPy's CSE functionality.
    
    No equation is given for this in IAPWS-95, and the derivative was 
    symbolically computed with SymPy.
    
    Like many higher-order derivatives of functions with exponentials, this one
    balloons to use many, many terms. Unfortunately, there are large numerical
    issues present at very low densities, under 0.0001 kg/m^3. No resolution
    to this issue besides switching to numerical derivatives or using arbitrary
    precision to compute it is known.
 
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-4
    kg/m^3 to 5000 kg/m^3, 90000 points were evaluated with mpmath. The mean 
    relative error was 9.767e-13, with a maximum relative error of 2.68E-10 and a 
    standard deviation of 7.38E-12.    
    
    Evaluating 40000 points in the 1e-10 to 1e-4 range, the 
    mean relative error was 1.14E-6, maximum relative error 0.000152, and
    standard deviation 6.59E-6. Caution is recommended in that range.
    
    Examples
    --------
    >>> iapws95_d3Ar_ddelta3(647.096/300.0, 999.0/322)
    0.33621190578672067
    '''
    ans = 0.0
    # Last two terms - appears to be working fine with mpmath
    x0 = delta - 1.0
    x1 = x0**2
    x2 = (tau - 1.0)**2
    x3 = x1**3.5
    x4 = delta - 1.0
    x5 = x4**2
    x6 = x5**1.66666666666666674
    x7 = -tau + 0.320000000000000007*x6 + 1.0
    x8 = 0.200000000000000011*x3 + x7**2
    x9 = 56.0*x1
    x10 = delta*x0
    x11 = 1/x8
    x12 = x3/x0
    x13 = x6/x5
    x14 = x7*(0.640000000000000013*delta - 0.640000000000000013)
    x15 = x13*x14
    x16 = x11*(1.40000000000000013*x12 + 3.33333333333333348*x15)
    x17 = delta*x16
    x18 = x3/x1
    x19 = x5**3.33333333333333348
    x20 = x4**(-4)
    x21 = x19*x20
    x22 = x1*x21
    x23 = x13*(-0.640000000000000013*tau + 0.20480000000000001*x6 + 0.640000000000000013)
    x24 = x20*x6
    x25 = x24*x7
    x26 = x1*x25
    x27 = x6/x4**3
    x28 = 1.28000000000000003*delta - 1.28000000000000003
    x29 = x28*x7
    x30 = x27*x29
    x31 = 0.42000000000000004*x12 + x15
    x32 = x11*x31**2
    x33 = 3.33333333333333348*x27
    x34 = 8.40000000000000036*x18 + 2.27555555555555555*x22 + 3.33333333333333348*x23 + 7.11111111111111072*x26 - x29*x33
    x35 = x10*x11
    x36 = x31**3/x8**2
    x37 = x11*(3.33333333333333348*x13*(-0.320000000000000007*tau + 0.102400000000000005*x6 + 0.320000000000000007) - x14*x33 + 4.20000000000000018*x18 + 1.13777777777777778*x22 + 3.55555555555555536*x26)
    x38 = 11.1111111111111107*x21
    x39 = x4**(-5)
    x40 = x1*x19
    x41 = (2.0*delta - 2.0)/x4**6
    x42 = x1*x6*x7
    x43 = 11.1111111111111107*x14*x24 + 11.1111111111111107*x25*x28 + 3.33333333333333348*x25*(3.83999999999999986*delta - 3.83999999999999986) - x33*(-2.56000000000000005*tau + 0.819200000000000039*x6 + 2.56000000000000005) - x37*(8.40000000000000036*x12 + 20.0*x15) + x38*(0.20480000000000001*delta - 0.20480000000000001) + x38*(0.40960000000000002*delta - 0.40960000000000002) - 13.6533333333333342*x39*x40 - 42.6666666666666643*x39*x42 + 11.3777777777777818*x40*x41 + 11.8518518518518547*x41*x42 + 42.0*x3/x0**3
    x44 = delta*x11
    x45 = 64.0*x1
    ans += (-0.148746408567240002*x8**0.849999999999999978*(9408.0*x1 - 3136.0*x10*(x9 - 3.0) + x11*(21.4199999999999982*x18 + 5.8026666666666662*x22 + 8.5*x23 + 18.1333333333333293*x26 - 8.5*x30 - 4.25000000000000266*x32) - x16*(285.599999999999966*delta - 285.599999999999966) + 142.799999999999983*x17*(x9 - 1.0) - 142.799999999999983*x35*(-1.66666666666666785*x32 + x34) + 0.849999999999999978*x44*(6.38888888888889994*x36 + x37*(7.14000000000000057*x12 + 17.0*x15) + x43) - 168.0)*exp(-28*x1 - 700*x2) + 0.318061108784439994*x8**0.949999999999999956*(12288.0*x1 - 4096.0*x10*(x45 - 3.0) + x11*(23.9399999999999977*x18 + 6.48533333333333228*x22 + 9.5*x23 + 20.2666666666666622*x26 - 9.5*x30 - 1.58333333333333259*x32) - x16*(364.799999999999955*delta - 364.799999999999955) + 182.399999999999977*x17*(x45 - 1.0) - 182.399999999999977*x35*(-0.555555555555555358*x32 + x34) + 0.949999999999999956*x44*(1.94444444444444287*x36 + x37*(7.97999999999999954*x12 + 19.0*x15) + x43) - 192.0)*exp(-32*x1 - 800*x2))


#    0-7 and 51-54 - works fairly fine. Checked with mpmath 
    # 1.5623946580944904e-16 6.3872154394719985e-15 2.1053159215966843e-12
    x0 = delta - 1.0
    x1 = x0**2
    x2 = -20.0*x1
    x3 = exp(x2 - 219.615*(0.826446280991736*tau - 1)**2)
    x4 = delta**3
    x5 = delta**2
    x6 = 6.0*delta - 240.0*x0*x5 + 1600.0*x1*x4 - 40.0*x4
    x7 = x3*x6
    x8 = delta*x0
    x9 = x1*x5
    x10 = 2.0*delta - 2.0
    x11 = x10*x4
    x12 = 1600.0*x11 - 360.0*x5 - 480.0*x8 + 4800.0*x9 + 6.0
    x13 = tau**4*exp(x2 - 390.625*(0.8*tau - 1)**2)
    ans += (0.211414783445121601*delta*tau - 0.046919851012788602*tau**0.375 - 630.92280475561995*tau*x10*x7 + 31.5461402377809996*tau*x12*x3 + 50426.3086833900015*x10*x13*x6 - 2521.31543416949989*x12*x13 + x3*(-50090.0165174960057*x11 + 11270.2537164366004*x5 + 15027.0049552487999*x8 - 150270.049552487995*x9 - 187.837561940610016) + x7*(1252.25041293740014*delta - 1252.25041293740014))

    x0 = 1/delta
    x1 = 1.33713144615930002 - 0.668565723079650009*delta
    x2 = tau**4
    x3 = -delta
    x4 = exp(x3)
    x5 = x2*x4
    x6 = tau**6
    x7 = 0.408676219019300013 - 0.204338109509650007*delta
    x8 = 0.000132425210079373988 - 0.0000662126050396869941*delta
    x9 = (x3 + 1.0)*(x3 + 2.0)
    x10 = tau*x4
    x11 = 0.160748684862510011*delta
    x12 = delta - 3.0
    x13 = x12*(delta - 2.0)
    x14 = tau**2
    x15 = 0.0801856578516139951*delta
    x16 = delta - 4.0
    x17 = x3 + 4.0
    x18 = x17*(x3 + 3.0)
    x19 = delta*x4
    x20 = tau**13
    x21 = 7.8686845206508003e-7*delta
    x22 = x12*x16
    x23 = delta**2
    x24 = tau**9
    x25 = x17*(x3 + 5.0)
    x26 = delta**4
    x27 = tau**3
    x28 = delta - 7.0
    x29 = x28*(delta - 6.0)
    x30 = delta**6
    x31 = delta - 9.0
    x32 = (x3 + 8.0)*(x3 + 9.0)
    x33 = delta**7
    x34 = delta - 10.0
    x35 = x31*x34
    x36 = delta - 11.0
    x37 = x34*x36
    x38 = (x3 + 12.0)*(x3 + 13.0)
    x39 = (x3 + 14.0)*(x3 + 15.0)
    x40 = 2.0*x23
    x41 = -x40
    x42 = x41 + 3.0
    x43 = -x23
    x44 = exp(x43)
    x45 = tau**7*x44
    x46 = 6.0*x23
    x47 = x43 + 1.0
    x48 = x47*(x41 + 1.0)
    x49 = delta*(x46 - x48 - 5.0)
    x50 = x24*x44
    x51 = tau**10
    x52 = x44*x51
    x53 = delta*x52
    x54 = 2.32333599943035995*x23
    x55 = -4.0*x23
    x56 = x42*x47
    x57 = x23*(x55 + 9.0)
    x58 = x43 + 2.0
    x59 = x42*x58
    x60 = x40 - (x23 - 2.0)*(x40 - 3.0)
    x61 = -x23*x60 + x57 + x60
    x62 = 1.91229319661920005*x23
    x63 = x41 + 5.0
    x64 = x58*x63
    x65 = delta**3
    x66 = x55 + 13.0
    x67 = 0.0821117635837920051*x23
    x68 = x43 + 3.0
    x69 = x63*x68
    x70 = x44*x6
    x71 = 0.545457404413720037*x23
    x72 = (x23 - 3.0)*(x40 - 5.0)
    x73 = 0.0567225376024680025*x23
    x74 = x68*(x41 + 7.0)
    x75 = x55 + 19.0
    x76 = 0.0333306019522852034*x23
    x77 = x41 + 9.0
    x78 = x77*(x43 + 4.0)
    x79 = x40 - x78
    x80 = x30*x44
    x81 = 0.116209344038340004*x23
    x82 = (x23 - 4.0)*(x40 - 9.0)
    x83 = x40 - x82
    x84 = 0.15446034229682401*x23
    x85 = 0.0815739460548159934*x23
    x86 = tau**8
    x87 = 0.00662162002549360022*x23
    x88 = x55 + 21.0
    x89 = 0.00798222879181639984*x23
    x90 = x77*(x43 + 5.0)
    x91 = x40 - x90
    x92 = 0.000634812332966279988*x23
    x93 = 0.0000655542733701200008*x23
    x94 = (x23 - 6.0)*(x40 - 11.0)
    x95 = -x65
    x96 = exp(x95)
    x97 = 0.392522541514298995*x65
    x98 = 3.0*x65
    x99 = -x98
    x100 = x99 + 4.0
    x101 = x95 + 1.0
    x102 = x101*(x99 + 2.0)
    x103 = 0.629892098347770046*x65
    x104 = x99 + 5.0
    x105 = x100*x101
    x106 = delta*x96
    x107 = tau**23
    x108 = (x65 - 1.0)*(x98 - 4.0)
    x109 = x100*x104
    x110 = 8.0*x26
    x111 = (2.0*x26 - 7.0)*(4.0*x26 - 13.0)
    x112 = 2.00560026836322003e-8*x30
    x113 = -12.0*x30
    x114 = (2.0*x30 - 1.0)*(3.0*x30 - 1.0)
    x115 = 6.0*x30
    x116 = -x30
    x117 = exp(x116)
    x118 = tau**50*x117
    x119 = x113 + 17.0
    x120 = 7.16605860758687996*x30
    x121 = (x115 - 5.0)*(x30 - 1.0)
    x122 = x115 - x121
    x123 = x117*x65
    x124 = 11.439899039065681*x30
    x125 = (5.0 - x115)*(x116 + 1.0)
    x126 = 4.26282567335316021*x30
    ans += (-x5*(x0*x1 + x0*(1.33713144615930002*delta - 1.33713144615930002) + x1))
    ans += (x4*x6*(x0*x7 + x0*(0.408676219019300013*delta - 0.408676219019300013) + x7))
    ans += (-tau**12*x4*(x0*x8 + x0*(0.000132425210079373988*delta - 0.000132425210079373988) + x8))
    ans += (-x10*(0.576981634680060029*delta - 0.192327211560020001*x9 - 0.769308846240080002))
    ans += (-tau**5*x4*(0.771271290103139995*delta - 0.25709043003437998*x9 - 1.02836172013751992))
    ans += (x5*(0.321497369725020021*delta*x12 + x11*(delta - x13) - x11 + 0.160748684862510011*x13))
    ans += (x14*x19*(-0.0400928289258069975*delta*(delta - x18) - x15*x16 + x15 - 0.0801856578516139951*x18))
    ans += (x19*x20*(3.93434226032540015e-7*delta*(delta - x22) + x16*x21 - x21 + 7.8686845206508003e-7*x22))
    ans += (x23*x24*x4*(-0.0000151882754176288*delta*(delta - 5.0) - 7.59413770881439999e-6*delta*(delta - x25) + 0.0000227824131264432008*delta - 0.0000227824131264432008*x25))
    ans += (x26*x27*x4*(0.00112501958703776009*delta*x28 + 0.000562509793518880044*delta*(delta - x29) - 0.00281254896759440022*delta + 0.00281254896759440022*x29))
    ans += (x30*x5*(-0.000031217304514269997*delta*x31 - 0.0000156086522571349985*delta*(delta - x32) + 0.00010926056579994499*delta - 0.00010926056579994499*x32))
    ans += (tau**11*x33*x4*(2.30759928459020003e-9*delta*x34 + 1.15379964229510002e-9*delta*(delta - x35) - 9.23039713836080014e-9*delta + 9.23039713836080014e-9*x35))
    ans += (delta**8*x5*(7.31643302884080017e-7*delta*x36 + 3.65821651442040008e-7*delta*(delta - x37) - 3.29239486297836013e-6*delta + 3.29239486297836013e-6*x37))
    ans += (delta**10*x20*x4*(-2.65023601493360004e-12*delta*(delta - 13.0) - 1.32511800746680002e-12*delta*(delta - x38) + 1.45762980821348004e-11*delta - 1.45762980821348004e-11*x38))
    ans += (delta**12*x10*(-1.25279173824907999e-9*delta*(delta - 15.0) - 6.26395869124539993e-10*delta*(delta - x39) + 8.14314629861901949e-9*delta - 8.14314629861901949e-9*x39))
    ans += (x45*(-0.431744036357279981*x23*x42 - 1.29523210907183994*x23 + 0.647616054535919972))
    ans += (0.0704459640350079963*tau*x44*x49)
    ans += (0.885291806701840045*x49*x50)
    ans += (1.60990679054111996*x53*(-x46 + x48 + 5.0))
    ans += (x52*(x54*(x40 - x56) - x54*(x55 + 7.0) - x54 + 1.16166799971517998*x56))
    ans += (0.0199876587963223988*delta*x27*x44*(x23*(x40 - x59) + x41 - x57 + x59))
    ans += (0.125434802850196009*delta*x45*x61)
    ans += (2.97263718841364*x53*x61)
    ans += (x23*x52*(-5.73687958985760016*x23 + x62*(x40 - x64) - x62*(x55 + 11.0) + 2.86843979492880008*x64))
    ans += (x65*x70*(-0.328447054335168021*x23 - x66*x67 + x67*(x40 - x69) + 0.16422352716758401*x69))
    ans += (x52*x65*(2.18182961765488015*x23 + x66*x71 - x71*(x40 - x72) - 1.09091480882744007*x72))
    ans += (x26*x52*(-0.28361268801234002*x23 + x73*(x40 - x74) - x73*(x55 + 15.0) + 0.14180634400617001*x74))
    ans += (tau*x80*(-0.233314213665996417*x23 - x75*x76 + x76*x79 + 0.116657106832998209*x78))
    ans += (x14*x80*(0.813465408268379986*x23 + x75*x81 - x81*x83 - 0.406732704134189993*x82))
    ans += (x27*x80*(-1.08122239607776804*x23 - x75*x84 + 0.540611198038884022*x78 + x79*x84))
    ans += (x2*x80*(0.571017622383712009*x23 + x75*x85 - 0.285508811191856005*x82 - x83*x85))
    ans += (x80*x86*(0.0463513401784551998*x23 + x75*x87 - 0.0231756700892275999*x82 - x83*x87))
    ans += (x33*x70*(-0.0638578303345311987*x23 - x88*x89 + x89*x91 + 0.0319289151672655994*x90))
    ans += (x33*x50*(-0.00507849866373023991*x23 - x88*x92 + 0.00253924933186511995*x90 + x91*x92))
    ans += (delta**9*x44*x86*(0.000655542733701200008*x23 - x93*(x40 - x94) + x93*(x55 + 25.0) - 0.000327771366850600004*x94))
    ans += (tau**16*x96*(-0.78504508302859799*x100*x65 + 0.130840847171432989*x102 + x97*(-x102 + x98) - x97))
    ans += (tau**22*x106*(-x103*x104 - x103 + 0.209964032782590015*x105 + 0.314946049173885023*x65*(-x105 + x98)))
    ans += (x106*x107*(1.38218756120317798*x104*x65 - 0.460729187067726031*x108 - 0.691093780601589103*x65*(-x108 + x98) + 1.38218756120317821*x65))
    ans += (x107*x23*x96*(0.067338831996018006*x109 + 0.067338831996018006*x65*(-x109 + 9.0*x65) - 1.21209897592832405*x65*(x95 + 2.0) - 0.606049487964162026*x65))
    ans += (delta**11*x51*(-0.00150455304995244002*x111 + 0.00100303536663496002*x26*(31.0 - x110) - 0.000501517683317480008*x26*(x110 - x111) + 0.0120364243996195202*x26)*exp(-x26))
    ans += (x118*(x112*(x113 + 11.0) - x112*(-x114 + x115) + x112 - 3.34266711393870005e-9*x114))
    ans += (tau**44*x123*(x119*x120 - x120*x122 - 4.77737240505792027*x121 + 28.6642344303475198*x30))
    ans += (tau**46*x123*(-x119*x124 + x124*(x115 - x125) + 7.62659935937712063*x125 - 45.7595961562627238*x30))
    ans += (x118*x65*(x119*x126 - 2.84188378223544014*x121 - x122*x126 + 17.0513026934126408*x30))
    return ans

def iapws95_dAr_dtau(tau, delta):
    r'''Calculates the first derivative of residual Helmholtz energy of water
    with respect to `tau` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    dAr_dtau : float
        Derivative of residual Helmholtz energy A/(RT) with respect to `tau`,
        [-]

    Notes
    -----
    This is an optimized implementatation taking 9 exp calls, 4 sqrts, 
    and 2 powers. It was generated using SymPy's CSE functionality, with a 
    limited amount of `horner` polynomial optimizations as well. It is over
    10x faster than a naive implementation.
    
    This implementation has been tested against a straightforward 
    implementation with the equations given in IAPWS-95.
    
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-10
    kg/m^3 to 5000 kg/m^3, 250000 points were evaluated. The mean relative
    error was 5.68E-14, with a maximum relative error of 6.73E-9 and a standard
    deviation of 1.35E-11.
    
    Over the same range, the model was evaluated to a precision of 50 
    decimal places with `mpmath`, and on 90000 points, the mean relative
    error was 4.66E-14, with a maximum relative error of 4.25E-10 and a standard
    deviation of 1.77E-12.
    
    The maximum error ocurs in the extremely low density regime,
    :math:`\rho < 1e-6`.

    Examples
    --------
    >>> iapws95_dAr_dtau(647.096/300.0, 999.0/322)
    -7.7043336309570
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    delta8 = delta4*delta4
    x1 = exp(-delta)
    x4 = exp(-delta2)
    x26 = exp(-delta3)
    x2 = delta*x1
    x5 = delta*x4
    x7 = x4*delta8
    x8 = 300.0*tau
    x9 = (delta - 1.0)
    x9 *= x9
    x10 = -20.0*x9
    x11 = (0.826446280991736*tau - 1.0)
    x11 = exp(x10 - 219.615*x11*x11)
    tau2 = tau*tau
    tau4 = tau2*tau2
    taurtinv = 1.0/sqrt(tau)
    tau_inv = taurtinv*taurtinv
    tau4inv = sqrt(taurtinv)
    tau8inv = sqrt(tau4inv)
    x15 = tau - 1.0
    x13 = tau*delta2
    x16 = x15*x15
    x17 = x9**1.66666666666666674
    x18 =  (-tau + 0.32*x17 + 1.0)
    x18 = 0.2*x9*x9*x9*sqrt(x9) + x18*x18
    x18_05 = x18**-0.05
    x19 = delta4*delta2
    x20 = delta3*x4
    x21 = delta8*delta2*x1
    x22 = delta4*delta*x4
    x24 = delta8*delta*x4
    x27 = delta*x26
    x50 = (0.8*tau - 1.0)
    return (delta*(-6.26395869124539993e-10*delta8*delta4*delta2*x1 + 0.159012546727090004*delta*taurtinv
             - 0.00626677396776150007*tau_inv*taurtinv - tau8inv*(0.196091503945184986*delta*tau8inv
            + 0.00293249068829928763*taurtinv*delta2 - 6.90879303824744984)
            - tau*(-tau*(tau*(tau*(tau*(tau*(tau*(-0.000131108546740240002*delta8*delta3*x4
            + tau*(tau*(-0.000626897104146849956*delta8*delta4*delta*exp(-delta4) + tau*(tau*(-0.000794551260476243984*x1
            + x13*(tau2*tau*(tau4*tau2*(tau*(tau4*tau4*tau4*tau4*tau4*tau*(tau2*(14.6176487721394821*delta3 - tau4*(5.92059121299049984*delta3
            + 2.78555592828224976e-8)) - 8.75851607593951975*delta3)*exp(-x19) + 0.516264378636137944*x26*delta2
            - 1.76612855042628292*x27) + 0.76986812020283002*x27) + 0.697817851580975979*x26) 
            + 5.11464493842302041e-6*x2 - 1.72265340970684017e-11*x21)) + 1.26917960652461008e-8*x1*delta8*delta)
            + 0.14180634400617001*x19*x4 - 7.43159297103409955*x20 - 1.36364351103430015*x22 
            + 4.78073299154799969*delta4*x4 + 5.80833999857589944*delta2*x4 - 4.02476697635280001*x5)
            - 0.0000683472393793296059*x1*delta4 + 0.00142832774917412992*x24 + 1.9919065650791401*x5)
            - 0.0132432400509872004*x7) - 0.219510904987843009*x20 - 0.755552063625239967*x4)
            + 1.22602865705789998*x1 + 0.123167645375688001*x22 + 0.0119733431877246006*x24)
            - 1.28545215017189984*x2) + 0.642994739450040043*x1*delta2 - 0.000062434609028539994*x1*delta8
            - 2.67426289231860004*x1 + 1.46328660576816003e-6*x21 - 0.0815739460548159934*x7) 
            + 0.00168752938055664013*x1*x19 + 0.0149907440972417982*x20 + 0.115845256722618001*x7)
            + 0.0801856578516139951*delta3*x1 + 0.058104672019170002*x7) + 0.00880894931021340005*delta3
            + 31.5461402377809996*x11*x13*(tau_inv - x8 + 363.0) + 31.3062603234350014*x11*delta2*(x8 - 363.0)
            - 2521.31543416949989*tau4*delta2*(-500.0*tau + 4.0*tau_inv + 625.0)*exp(x10 - 390.625*x50*x50)
            - 0.192327211560020001*x2 + 0.0176114910087519991*x5 + 0.00833265048807130086*x7 
            + 0.148746408567240002*(1400.0*x15*x18*x18_05*x18_05*x18_05
            + x18_05*x18_05*x18_05*(-1.7*tau + 0.544*x17
            + 1.7))*exp(-700.0*x16 - 28.0*x9) - 0.318061108784439994*(1600.0*x15*x18*x18_05 
            + x18_05*(-1.9*tau + 0.608*x17 + 1.9))
            *exp(-800.0*x16 - 32.0*x9) - 8.78032033035609949))

def iapws95_d2Ar_dtau2(tau, delta):
    r'''Calculates the second derivative of residual Helmholtz energy of water
    with respect to `tau` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    d2Ar_dtau2 : float
        Second derivative of residual Helmholtz energy A/(RT) with respect to
        `tau`, [-]

    Notes
    -----
    This is an optimized implementatation taking 9 exp calls, 4 sqrts, 
    and 2 powers. It was generated using SymPy's CSE functionality, with a 
    limited amount of `horner` polynomial optimizations as well. It is over
    10x faster than a naive implementation.
    
    This implementation has been tested against a straightforward 
    implementation with the equations given in IAPWS-95.
    
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-10
    kg/m^3 to 5000 kg/m^3, 4E6 points were evaluated. The mean relative
    error was 4.595E-16, with a maximum relative error of 1.835e-10 and a standard
    deviation of 1.209E-13.
    
    Over the same range, the model was evaluated to a precision of 50 
    decimal places with `mpmath`, and on 250000 points, the mean relative
    error was 2.6026E-16, with a maximum relative error of 2.36E-12 and a standard
    deviation of 8.055E-15.

    This comparison indicates this implementation is more accurate than the
    straightforward implementation.
    
    Examples
    --------
    >>> iapws95_d2Ar_dtau2(647.096/300.0, 999.0/322)
    -1.2616419775539
    '''
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta*delta3
    delta5 = delta*delta4
    delta6 = delta*delta5
    delta8 = delta2*delta6
    delta9 = delta8*delta
    x1 = exp(-delta)
    x3 = x1*delta3
    tau2 = tau*tau
    tau4 = tau2*tau2
    tau6 = tau4*tau2
    tau7 = tau*tau6
    tau8 = tau*tau7
    tau11 = tau8*tau2*tau
    tau21 = tau11*tau8*tau2
    tau48 = tau21*tau21*tau6
    x5 = x1*tau2
    x8 = exp(-delta2)
    taurtinv = 1.0/sqrt(tau)
    tau_inv = taurtinv*taurtinv
    tau4inv = sqrt(taurtinv)
    tau8inv = sqrt(tau4inv)
    x10 = x8*delta8
    x11 = tau*tau4*x8
    x13 = delta*x8
    x15 = delta3*x8
    x19 = tau8*x8
    x20 = exp(-delta3)
    x22 = exp(-delta6)
    x23 = delta3*x20
    x26 = x22*delta5
    x28 = delta9*x8
    x29 = 0.82644628099173556*tau
    x30 = (x29 - 1.0)
    x30 *= x30
    x31 = (delta - 1.0)
    x31 *= x31
    x32 = -20.0*x31
    x33 = delta2*exp(-219.615*x30 + x32)
    x34 = 1.0/tau2
    x35 = 1.0/tau
    x36 = 0.8*tau
    x37 = tau - 1.0
    x38 = x37*x37
    x39 = x31**1.66666666666666674
    x40 = (-tau + 0.32*x39 + 1.0)
    x40 *= x40
    x41 = 0.2*x31*x31*x31*sqrt(x31) + x40
    x43 = x41**(-0.05)
    x42 = x43*x43*x43
    
    x50 = 1.0 - x29 + 0.00275482093663911853*x35
    x51 = (0.0064*x35 - x36 + 1.0)
    return (delta*(delta9*delta*(-0.00564207393732165004*delta3*tau8*exp(-delta4)
            - 2.0671840916482082e-10*delta2*x1*tau11 - 0.000917759827181680011*delta*tau6*x8
            + 4.38985981730447989e-6*x5 )
            - 5.14180860068759937*delta*tau2*tau*x1 
            +tau_inv*(
            + 0.0490228759862962465*delta*tau4inv
            +taurtinv*(- 0.079506273363545002*delta
            + 0.00940016095164225053*tau_inv
            + 0.00183280668018705488*tau8inv*delta2)
            - 0.863599129780931229*tau8inv)
            
            + tau8*tau*(tau11*(tau11*tau11*x26*(657.79419474627673*tau2 
            - 376.616191265399323) + 16.1672305242594305*x23)
            + tau*(10.4672677737146405*tau4*delta2*x20 - 0.00874006386523868382*x1)
            + 1.26917960652460995e-7*x1*delta9)
            
            + tau*(0.00337505876111328026*x1*delta6
            + 0.231690513445236002*x10 + 0.0299814881944835965*x15 
            - 31.5461402377809996*x33*(x34 - 131769.0*x50*x50 
            + 300.0)) + 52.2750599871830985*delta2*x19 - 1.3649224048583023e-6*delta2*tau48*x22
            + 1.92898421835012002*delta2*x5 + 2521.31543416949989*delta2*tau4*(4.0*x34 - 390625.0*x51*x51
            + 500.0)*exp(x32 - 390.625*(x36 - 1.0)*(x36 - 1.0)) - 0.000546777915034636847*x1*tau7*delta4
            + 6.13014328528949992*x1*tau4 - 0.0927026803569103997*x10*tau6 - 0.24472183816444798*x10*tau2
            - 0.058104672019170002*x10 - 1.31706542992705811*x11*delta3 - 4.5333123817514398*x11
            + 15.9352525206331208*tau7*x13 + 0.0114266219933930394*tau7*x28 - 36.222902787175201*x13*tau8
            - 66.8843367393068888*tau8*x15 + 0.0000613757392610762517*tau11*x3 + 43.0265969239319972*delta4*x19
            + 11.3578163299950354*delta4*x20*tau21 - 12.272791599308702*x19*delta5 + 1.27625709605553017*x19*delta6
            - 290.108969436534494*tau48*x26 - 38.8548281093782251*x23*tau21 + 0.615838226878440032*delta5*tau4*x8
            + 0.0598667159386229997*x28*tau4 - 0.0801856578516139951*x3 - 31.3062603234350014*x33*(131769.0*x30 - 300.0)
            - 0.000187303827085619982*x5*delta8 - 8.02278867695580011*x5
            - 0.148746408567240002*(2800.0*x37*x42*(-1.7*tau 
            + 0.544*x39 + 1.7) - 0.51*x40*x42/(x41) 
            + 1400.0*x41*x42*(1400.0*x38 - 1.0) + 1.7*x42)*exp(-28.0*x31 - 700.0*x38) 
            + 0.318061108784439994*(3200.0*x37*x43*(-1.9*tau + 0.608*x39 
            + 1.9) - 0.19*x40*x43/(x41)
            + 1600.0*x41*x43*(1600.0*x38 - 1.0) + 1.9*x43)*exp(-32.0*x31 - 800.0*x38)))


def iapws95_d2Ar_ddeltadtau(tau, delta):
    r'''Calculates the second derivative of residual Helmholtz energy of water
    with respect to `tau` and also `delta` according to the IAPWS-95 standard.

    Parameters
    ----------
    tau : float
        Dimensionless temperature, (647.096 K)/T [-]
    delta : float
        Dimensionless density, rho/(322 kg/m^3), [-]

    Returns
    -------
    d2Ar_ddeltadtau : float
        Second derivative of residual Helmholtz energy A/(RT) with respect to
        `tau` and `delta`, [-]

    Notes
    -----
    This is an optimized implementatation taking 11 exp calls, 4 sqrts, 
    and 3 powers. It was generated using SymPy's CSE functionality, with 
    select polynomial optimizations by hand as well. It is over
    10x faster than a naive implementation.
    
    This implementation has been tested against a straightforward 
    implementation with the equations given in IAPWS-95.
    
    Over a linear temperature range of 200 K to 5000 K and a logarithmic 
    density range of 1E-10
    kg/m^3 to 5000 kg/m^3, 4E6 points were evaluated. The mean relative
    error was 2.82E-14, with a maximum relative error of 8.404E-9 and a standard
    deviation of 5.166e-12.    
    
    Over the same range, the model was evaluated to a precision of 50 
    decimal places with `mpmath`, and on 90000 points, the mean relative
    error was 6.974E-14, with a maximum relative error of 4.286E-9 and a standard
    deviation of 4.286E-12.
    
    Examples
    --------
    >>> iapws95_d2Ar_ddeltadtau(647.096/300.0, 999.0/322)
    -0.198403562385
    '''
    delta2 = delta*delta
    delta3 = delta2*delta
    delta4 = delta3*delta
    delta5 = delta4*delta
    delta6 = delta5*delta
    delta8 = delta6*delta2
    delta9 = delta8*delta
    tau2 = tau*tau
    tau3 = tau*tau2
    tau4 = tau*tau3
    tau5 = tau*tau4
    tau7 = tau2*tau5
    tau8 = tau*tau7
    tau9 = tau*tau8
    tau16 = tau8*tau8
    tau32 = tau16*tau16
    tau49 = tau32*tau*tau16

    taurtinv = 1.0/sqrt(tau)
    tau_inv = taurtinv*taurtinv
    tau4inv = sqrt(taurtinv)
    tau8inv = sqrt(tau4inv)

    x2 = exp(-delta)
    x10 = exp(-delta2)
    x24 = exp(-delta3)
    x29 = exp(-delta6)
    
    x3 = delta*x2*(delta - 2.0)
    x5 = delta - 1.0
    x6 = x2*x5
    x9 = delta3*(delta - 4.0)
    x11 = delta*x10*(delta2 - 1.0)
    x12 = x2*tau3
    x13 = tau7*tau5*x2
    x21 = 2.0*delta2
    x22 = x10*delta8*(x21 - 9.0)
    x23 = tau3*tau3*x10
    x25 = delta3*(delta2 - 2.0)
    x26 = x10*tau9
    x28 = delta5*(delta2 - 3.0)
    x30 = delta5*x29*(delta6 - 1.0)
    x32 = x10*delta9*(delta2 - 5.0)
    x34 = 3.0*delta3
    x35 = delta3*(x34 - 4.0)
    x36 = tau9*tau9*tau4*x24
    x37 = 300.0*tau
    x38 = x5*x5
    x38rt = sqrt(x38)
    x39 = -20.0*x38
    x40 = delta3*(-40.0*delta + 40.0 + 3.0/delta)
    x41 = (0.826446280991736*tau - 1.0)
    x41 = x40*exp(x39 - 219.615*x41*x41)
    x42 = 1.0/tau
    x57 = x38**0.666666666666666666
    x43 = x38*x57
    x44 = 1.0 -tau + 0.32*x43
    x44 = 0.2*x38*x38*x38*x38rt + x44*x44
    x45 = delta*x5
    x46 = 56.0*x45
    x47 = tau - 1.0
    x48 = x47*x47
    x49 = exp(-28.0*x38 - 700.0*x48)
    x50 = x47*x49
    x51 = 64.0*x45
    x52 = exp(-32.0*x38 - 800.0*x48)
    x53 = x44**(-0.05)
    x54 = x38 + 25.0*x48
    x55 = x53*x53*x53
    x56 = 1.4*x38*x38*x38rt
    x58 = -2.1333333333333333*tau + 0.682666666666666755*x43 + 2.1333333333333333
    x59 = x56 + x57*x58
    x60 = x45*x52
    x61 = x56 + x58*(x5*x5)**0.666666666666666741
    x100 = 0.8*tau - 1.0
    return (delta9*delta*(6.26395869124539993e-10*delta4*x2*(delta - 15.0) 
            + 0.00125379420829369991*delta3*tau9*(2.0*delta4 - 7.0)*exp(-delta4)
            + 1.72265340970684017e-11*delta2*x13*(delta - 13.0) 
            + 0.000262217093480480003*delta*x10*tau7*(delta2 - 6.0)
            - 1.46328660576816003e-6*x12*(delta - 11.0) )
                
            - 0.392183007890369972*delta*tau4inv
            + taurtinv*(0.318025093454180008*delta
            - 0.00626677396776150007*tau_inv
            - 0.00879747206489786202*tau8inv*delta2)
            + 6.90879303824744984*tau8inv
            
            + tau9*tau*(# These two terms do nto want to have tau32 factored out
                    -87.7058926328368926*tau3*tau32*x30 
            + 52.5510964556371221*tau*tau32*x30
            
            - 0.76986812020283002*tau9*tau2*x24*x35 
            - 2.09345355474292782*tau5*delta2*x24*(delta3 - 1.0)
            + 0.000794551260476243984*tau*x6 
            - 1.26917960652461008e-8*delta9*x2*(delta - 10.0))
            
            + 0.0801856578516139951*tau*x2*x9
            + 0.058104672019170002*tau*x22 + 31.5461402377809996*tau*x41*(-x37 + x42 + 363.0)
            + 0.0352357972408536002*delta3 - 0.642994739450040043*delta2*x12*(delta - 3.0)
            - 5.80833999857589944*delta2*x26*(x21 - 3.0)
            + 8.35666778484674929e-8*delta2*x29*tau49*(2.0*delta6 - 1.0) - 0.0299814881944835965*x10*tau2*x25 
            - 0.246335290751376002*x10*x28*tau5 - 3.9838131301582802*x11*tau8 + 8.04953395270560002*x11*tau9
            - 0.0352229820175039982*x11 + 0.000062434609028539994*x12*delta8*(delta - 9.0) 
            - 5.11464493842302041e-6*x13*x9 + 0.0000683472393793296059*delta4*tau8*x2*(delta - 5.0)
            - 4.78073299154799969*delta4*x26*(x21 - 5.0) - 0.516264378636137944*delta4*x36*(x34 - 5.0)
            - 0.00285665549834825984*tau8*x32 - 0.00168752938055664013*delta6*tau2*x2*(delta - 7.0)
            - 0.14180634400617001*delta6*x26*(x21 - 7.0) - 0.115845256722618001*tau2*x22
            + 0.0132432400509872004*x22*tau7 + 0.0815739460548159934*x22*tau3 - 0.00833265048807130086*x22 
            + 0.439021809975686017*x23*x25 + 0.755552063625239967*x23*(x21 - 1.0) + 14.8631859420681991*x25*x26 
            + 2.7272870220686003*x26*x28 + 1.28545215017189984*x3*tau4 + 0.192327211560020001*x3 
            + 35.5235472779429955*x30*tau49 - 0.0239466863754492013*x32*tau5 + 1.76612855042628292*x35*x36 
            + 2.67426289231860004*tau3*x6 - 2521.31543416949989*x40*tau4*(-500.0*tau + 4.0*x42 + 625.0)*exp(
            x39 - 390.625*x100*x100) + 31.3062603234350014*x41*(x37 - 363.0) 
            - 208.244971994135994*x44*x53*x53*x53*x50*(x46 - 1.0) 
            + 508.89777405510398*x44*x53*x47*x52*(x51 - 1.0) 
            - 0.148746408567240002*x45*x49*(x53*x53*x53/x44*x61*(
            -0.255*tau + 0.0816*x43 + 0.255) 
            - 1.81333333333333346*x55*x57) + 177.008226195015595*x45*x50*x55*x59 
            - 483.452885352348801*x47*x53*x59*x60 - 0.318061108784439994*x53*(-x51*x52 
            + exp(-32.0*x54))*(-1.9*tau + 0.608*x43 + 1.9)
            + 0.148746408567240002*x55*(-x46*x49 + exp(-28.0*x54))*(-1.7*tau
            + 0.544*x43 + 1.7) - 1.22602865705789998*x6*tau5
            + 0.318061108784439994*x60*(x53/x44*x61*(-0.095*tau 
            + 0.0304*x43 + 0.095) - 2.02666666666666666*x53*x57)
            - 8.78032033035609949)

### Vapor pressure solution
def _P_G_dG_dV_T_dG_dV_T(T, V):
    '''For calculating vapor pressure'''
    _MW_kg = MW/1000
    R_MW_0_001 = _MW_kg*R95
    rho = rho_mass = MW/V/1000
    tau = Tc/T
    delta = rho_mass/rhoc

    tau = Tc/T
    delta = rho/rhoc
    
    A = iapws95_A0(tau, delta) + iapws95_Ar(tau, delta)
    dA_dtau = iapws95_dAr_dtau(tau, delta) + iapws95_dA0_dtau(tau, delta)
    d2A_dtau2 = iapws95_d2Ar_dtau2(tau, delta) + iapws95_d2A0_dtau2(tau, delta)
    dA_ddelta = iapws95_dAr_ddelta(tau, delta) + 1/delta
    d2A_ddelta2 =  iapws95_d2Ar_ddelta2(tau, delta) -1/(delta*delta)
    d2A_ddeltadtau = iapws95_d2Ar_ddeltadtau(tau, delta)
    
    P = (dA_ddelta*delta)*rho*R95*T

    S = (tau*dA_dtau - A)*R_MW_0_001
    H = (tau*dA_dtau + dA_ddelta*delta)*T*R_MW_0_001
    G = H - T*S
    
    dP_dV = (-MW*MW*R95*T*(MW*d2A_ddelta2/10**9 + V*rhoc*dA_ddelta/500000)/(V*V*V*V*rhoc*rhoc))
    
    dS_dV_T = ((dA_ddelta - Tc*d2A_ddeltadtau/T)
                                   *R_MW_0_001*MW/1000/(V*V*rhoc))
    
    dS_dP_T = dS_dV_T/dP_dV
    
    dH_dV_T = (T*MW*(-MW*d2A_ddelta2/(V*rhoc*10**6) - dA_ddelta/1000
         -   Tc*d2A_ddeltadtau/T/1000)*R_MW_0_001/(rhoc*V*V))
    dH_dP_T = dH_dV_T/dP_dV
    
    dG_dP_T = -T*dS_dP_T + dH_dP_T
    
    dG_dV_T = dG_dP_T*dP_dV
    
    return G, P, dG_dV_T, dP_dV

def iapws95_sat_err_and_jac(Vs, T):
    #print(Vs, 'Vs')
    V_l, V_g = Vs
    G_l, P_l, dG_dV_l, dP_dV_l =  _P_G_dG_dV_T_dG_dV_T(T, V_l)
    G_g, P_g, dG_dV_g, dP_dV_g =  _P_G_dG_dV_T_dG_dV_T(T, V_g)
    err = [0.0]*2
    err[0] = G_l - G_g
    err[1] = P_l - P_g
    
    jac = [[dG_dV_l, -dG_dV_g],
           [dP_dV_l, -dP_dV_g]]
    #print(err, 'err')
#     print([float(i) for i in err])
    return err, jac

def iapws95_saturation(T, xtol=1e-16):
    # rhol_sat_IAPWS95
    rhog, rhol = rhog_sat_IAPWS(T), rhol_sat_IAPWS(T)
    Vg = MW/(rhog*1000)
    Vl = MW/(rhol*1000)
    (Vg, Vl), iters = newton_system(iapws95_sat_err_and_jac, [Vg, Vl], args=(T,),
                                    jac=True, xtol=xtol, solve_func=solve_2_direct,
                                    damping_func=damping_maintain_sign
                                    )
    rhol = MW/(Vl*1000)
    rhog = MW/(Vg*1000)
    Psat = iapws95_P(T, rhol)
    return Psat, rhol, rhog


Psat_coeffs_iapws95_235_273 = [-0.0404815673828125, 0.0458526611328125, 0.375885009765625, -0.4257164001464844, -1.6131267547607422, 1.8266944885253906, 4.242580890655518, -4.803114175796509, -7.645859479904175, 8.652797102928162, 10.002255886793137, -11.31277510523796, -9.817383006215096, 11.092894461005926, 7.370591592043638, -8.314844283275306, -4.277057046769187, 4.812127415090799, 1.925885249627754, -2.1572029151720926, -0.6719730040349532, 0.7471348317922093, 0.18054040489369072, -0.19828112466348102, -0.03694033686770126, 0.03975118442394887, 0.005655801688135398, -0.005886636653485766, -0.000631581975426343, 0.0006184879482020733, 5.696264696553044e-05, -4.779557498579834e-05, -1.56374487509936e-05, 5.210719629289429e-05, -0.0001258797167196235, 0.00031015444493132094, -0.0012805056242366497, 0.01149752302139273, -0.13780084439128548, 1.6313925771945534, -11.998047252356896]
Psat_coeffs_iapws95_273_460 = [8.147908374667168e-07, -1.3258541002869606e-06, -2.960063284263015e-06, 4.862929927185178e-06, 6.120302714407444e-06, -9.979732567444444e-06, -6.454627509810962e-06, 1.1027086657122709e-05, 5.648329533869401e-06, -1.0190725788561394e-05, -1.3649275842908537e-06, 4.478973778532236e-07, 2.4834392547745665e-05, -5.852100520087333e-05, 0.00011126565522801002, -0.00036335083247251987, 0.0013489128931727379, -0.0047478345366922525, 0.018241801868318586, -0.07380538465894615, 0.2868986774488336, -1.015108797564658, 3.4736507508968284, -5.619649447193197]
Psat_coeffs_iapws95_460_609 = [2.0064180716872215e-07, 2.892629709094763e-07, -8.132046787068248e-07, -1.2402815627865493e-06, 1.2566961231641471e-06, 1.8773480405798182e-06, -1.3796370694763027e-06, -1.8621731214807369e-06, 4.6375271267606877e-07, 1.1436318345658947e-06, 1.6705293433005863e-06, 1.7150333064819279e-06, 1.3039190918107124e-06, 5.15597523786937e-06, 1.185603442621641e-05, 1.694959209075364e-05, 5.210011446887819e-05, 7.716156925408058e-06, 0.0009684942799086382, -0.0037996536704614225, 0.031002148395983753, -0.1716204979856538, 1.2065145338661658, -1.519520204731372]
Psat_coeffs_iapws95_609_643 = [1.0573770850896835e-05, 1.1763535439968109e-05, -6.56498596072197e-05, -7.157295476645231e-05, 0.0001856974558904767, 0.00019913158030249178, -0.0003122607449768111, -0.00033338970388285816, 0.0003383217645023251, 0.00036499839916359633, -0.00024304199087055167, -0.0002661161202013318, 0.00012133384302615013, 0.00013108684390772396, -4.926522541381928e-05, -5.498319717389677e-05, 3.840095160967394e-06, 9.275973251732239e-06, 4.700250627820424e-06, 1.7356479738772634e-05, 3.5141797559923305e-05, 4.994296202665005e-05, 5.6967466662061206e-05, 5.500573225036831e-05, 7.36448261119127e-05, 0.000486011140590506, -0.0032944327011136323, 0.20513406870743833, -0.24748976782381205]
Psat_coeffs_iapws95_643_646 = [9.057896477315808e-08, 1.211502649312024e-07, -3.360720484124613e-07, -4.495756513733795e-07, 5.768298478869838e-07, 7.768647307671017e-07, -5.138533794024625e-07, -6.86156649720715e-07, 3.3836656099239804e-07, 4.158511894836181e-07, -1.3213676086643034e-07, -1.300110328283921e-07, 2.7710538380576466e-07, 8.91687420706555e-07, 1.834557936633563e-06, 2.582291823115257e-06, 1.5351452190368736e-06, 8.540061517871983e-07, 1.7006021977569903e-05, 0.02017998616540387, -0.02472793196736738]
Psat_coeffs_iapws95_near_critical = [-0.007006160914897919, -0.006299436092376709, 0.06787101551890373, 0.06092516239732504, -0.3032621801830828, -0.27187402336858213, 0.8286682133330032, 0.7421519589261152, -1.5481375773670152, -1.3854497998690931, 2.0945177564717596, 1.8733694699149055, -2.1206636809074553, -1.8960203259257469, 1.6377062088922685, 1.4638579778784333, -0.9743636551480677, -0.8708100283029125, 0.44806853431398963, 0.40042604495646117, -0.15886602644877712, -0.14197406009749258, 0.04308180067010703, 0.0385020246406107, -0.008810816846762659, -0.007874311835966985, 0.0013301372745293527, 0.0011886957505600204, -0.00014366835734647143, -0.0001283599871508001, 1.0624248652432722e-05, 9.493371248271015e-06, -4.747152791886172e-07, -4.074279138207085e-07, 7.613360554058896e-08, 1.2227703013014193e-07, 2.2590783514621522e-07, 6.263045703604046e-07, 2.5993277169989305e-06, 0.0022606230644268894, -0.00226439777925055]

def iapws95_Psat(T):
    # Fit to under 1e-12 precision, with the EOS evaluated with mpmath for max
    # precision.
    if 235.0 <= T < 273.15:
        # Equation solves down to this temperature but not below.
        coeffs = Psat_coeffs_iapws95_235_273
        val = horner(coeffs, 0.052424639580602915*(T - 254.074999999999989))
    elif 273.15 <= T <= 460.1225: # half the points
        coeffs = Psat_coeffs_iapws95_273_460
        val = horner(coeffs, 0.0106967602187487444*(T - 366.636250000000018))
    elif 460.1225 < T <= 609.7005:
        coeffs = Psat_coeffs_iapws95_460_609
        val = horner(coeffs, 0.0133709502734359297*(T - 534.911500000000046))
    elif 609.7005 < T <= 643.35555:
        coeffs = Psat_coeffs_iapws95_609_643
        val = horner(coeffs, 0.0594264456597153254*(T - 626.528025000000071))
    elif 643.35555 < T <= 646.721055:
        coeffs = Psat_coeffs_iapws95_643_646
        val = horner(coeffs, 0.594264456597155322*(T - 645.038302499999986))
    elif T < Tc:
        coeffs = Psat_coeffs_iapws95_near_critical
        val = horner(coeffs, 5.33411567172122325*(T - 646.908527499949969))
    else:
        raise ValueError("Temperature range must be between 273.15 K to 647.096 K")
    return exp(val)*22064000.0

def rhol_sat_IAPWS95(T):
    if 235.0 <= T < 273.15:
        coeffs = [-17.222183227539062, 18.77007293701172, 160.3718719482422, -174.77803325653076, -689.9691653251648, 751.9200625419617, 1818.6332448720932, -1981.8576150536537, -3283.683091402054, 3578.2904051095247, 4302.269122205675, -4688.138357363641, -4227.204455545172, 4606.233353788964, 3174.7852028068155, -3459.378792709904, -1840.9434151535388, 2005.933489420975, 826.9075718919048, -901.0039060153795, -287.0057090885457, 312.71981305931695, 76.35687688595681, -83.19737945724955, -15.353652797617087, 16.729144014339, 2.283299787289934, -2.4879351324522645, -0.2437405315180854, 0.2654055681062317, 0.01792607828713244, -0.019693974261834057, -0.0004738133771589048, 0.0007854872751575925, -0.0006095267475783039, 0.0022049570385940243, -0.005221612771147854, 0.01096251463881881, -0.026019484252807534, 0.036332822379610566, 3.0872453741478227]
        val = horner(coeffs, 0.052424639580602915*(T - 254.074999999999989))
    elif 273.15 <= T <= 460.1225: # half the points
        coeffs = [1.7102574929594994e-06, -3.5461271181702614e-06, -4.744477337226272e-06, 1.1568277841433883e-05, 1.0390001989435405e-05, -2.5427980290260166e-05, -7.063001248752698e-06, 2.3168839106801897e-05, 1.3829438557877438e-05, -3.149418807879556e-05, 1.2095285001123557e-05, -1.292190836466034e-05, 3.8989083236629085e-05, -6.452956441194146e-05, 9.861656508292072e-05, -0.00016529921067132136, 0.00028210334725997654, -0.0004912796557263732, 0.0008216553207591737, -0.0015053379092586638, 0.002406131617918583, -0.006780978799273418, 0.00908602594184732, -0.06335671046246652, -0.19965126834113162, 2.9904548603014143]
        val = horner(coeffs, 0.0106967602187487444*(T - 366.636250000000018))
    elif 460.1225 < T <= 609.7005:
        coeffs = [2.034008502960205e-06, 2.436339855194092e-06, -1.459755003452301e-05, -1.7838552594184875e-05, 4.556635394692421e-05, 5.651172250509262e-05, -8.218304719775915e-05, -0.00010129279689863324, 9.901652811095119e-05, 0.00011561789142433554, -9.266487904824317e-05, -0.00010660312182153575, 4.007856477983296e-05, 3.6490591810434125e-05, -5.036522770751617e-05, -5.477863123815041e-05, -3.106556266629923e-05, -4.485767993855916e-05, -7.495337143836878e-05, -0.00011451218220770443, -0.0001927638574699131, -0.00035594329908761324, -0.0006906136151680897, -0.0013645509703152925, -0.0026981001401753524, -0.005280856348563123, -0.01113746163068241, -0.024864014260847, -0.07935231088522635, -0.3683075956298905, 2.424989140035793]
        val = horner(coeffs, 0.0133709502734359297*(T - 534.911500000000046))
    elif 609.7005 < T <= 643.35555:
        coeffs = [-0.005523681640625, -0.00545501708984375, 0.0537109375, 0.05365753173828125, -0.2407855987548828, -0.24350690841674805, 0.6603143215179443, 0.6763124465942383, -1.2398340702056885, -1.2859440445899963, 1.6920538246631622, 1.7755838260054588, -1.7406707555055618, -1.845441060140729, 1.3810706129297614, 1.4778835410252213, -0.8551894403062761, -0.9241933110170066, 0.4143857065355405, 0.4525773521454539, -0.158058475019061, -0.17400417243334232, 0.047671283817180665, 0.05314236543927109, -0.010313240202776797, -0.01200570931268885, 0.0015920358180210314, 0.0011764186526761478, -0.0018841396051811898, -0.0023740562720035996, -0.002218709652339612, -0.001974439023143759, -0.0013470066813248138, -0.00040877893032842394, 0.00047615294604597125, 0.0007077648942215298, -0.0006431129948183401, -0.004859799769289075, -0.013341369671968695, -0.027982763294779223, -0.059743791257592174, -0.2308678349483143, 1.7406588358235677]
        val = horner(coeffs, 0.0594264456597153254*(T - 626.528025000000071))
    elif 643.35555 < T <= 646.721055:
        coeffs = [-0.0025634765625, -0.0032825469970703125, 0.022905349731445312, 0.029368877410888672, -0.09517884254455566, -0.12216567993164062, 0.24343305826187134, 0.3127277195453644, -0.42805875837802887, -0.5503446385264397, 0.5477428995072842, 0.7048549354076385, -0.5267451787367463, -0.678729840554297, 0.3876991346478462, 0.5006899123545736, -0.2203993159928359, -0.28577900753589347, 0.09690586793294642, 0.12655803375673713, -0.032833380970259896, -0.043427356979009346, 0.00839818057283992, 0.01135846923398276, -0.0016914226212065842, -0.0023644802264470854, 0.00010811350080075499, 0.00018391106829085402, -0.00024146455655227328, -0.00037838730947470367, -0.0005190572816360373, -0.0007225867524383034, -0.0009020460988393708, -0.0010130621772712622, -0.0011647964785088671, -0.0018405530481413468, -0.004267423867078957, -0.010988133931514904, -0.027298278157160433, -0.0959238698681941, 1.3178562243402536]
        val = horner(coeffs, 0.594264456597155322*(T - 645.038302499999986))
    elif 646.721055 < T <= 647.07:
        coeffs = [-0.0229644775390625, -0.03644561767578125, 0.212982177734375, 0.345855712890625, -0.9255542755126953, -1.5350914001464844, 2.500049591064453, 4.224307537078857, -4.699752390384674, -8.060362190008163, 6.527629733085632, 11.303938150405884, -6.944240380078554, -12.053825380280614, 5.790451696142554, 9.975626651197672, -3.842256080592051, -6.484273410635069, 2.048077084153192, 3.329579567653127, -0.8812978177811601, -1.3517990967593505, 0.3061402047842421, 0.432271851645055, -0.08556108074390067, -0.10814184474293143, 0.01886623827419953, 0.02070108844236529, -0.0034511384185407223, -0.003249207452029168, 0.00016951103771134512, -3.8158334660920445e-06, -0.00046348517059513483, -0.000568669470194294, -0.0007047358201628384, -0.000932012451160702, -0.0012536221138755854, -0.0017243174367475023, -0.0024525617113163867, -0.0036673798893029352, -0.005942076328912957, -0.011318687131713495, -0.03982414332624489, 1.1397375425430707]
        val = horner(coeffs, 5.73156228058745576*(T - 646.895527500000071))
    elif 647.07 < T <= 647.09:
        coeffs = [-1.2898817658424377e-07, -2.3213215172290802e-07, 4.318426363170147e-07, 1.466367393732071e-06, 1.237931428477168e-07, -4.1386374505236745e-06, -3.063873009523377e-06, 6.612863217014819e-06, 7.028365871519782e-06, -6.654002390860114e-06, -8.305772780659026e-06, 4.523619054452865e-06, 5.796317054773681e-06, -3.283294972789008e-06, -5.114287745300317e-06, -2.4331531847110455e-06, -3.3224186708480374e-06, -5.743677981229212e-06, -7.499835216862039e-06, -1.1263549756890257e-05, -2.327157668147173e-05, -5.8298991014504864e-05, -0.00015448293805092606, -0.00041322408313870795, -0.0011335557976513144, -0.003412156416544368, -0.014736997602476914, 1.0571055050155371]
        val = horner(coeffs, 100.000000000090949*(T - 647.080000000000041))
    elif 647.09 < T <= 647.095:
        # On 500 points, avg err 1.0048621666898271e-13, stdev 9.80697435140731e-14, max err 5.999286627611623e-13
        # Simply can't get better
        coeffs = [-0.0005059242248535156, -0.00014829635620117188, 0.004393354058265686, 0.0012152865529060364, -0.017519604414701462, -0.004545917734503746, 0.04250476974993944, 0.010283735115081072, -0.07004969473928213, -0.015714182169176638, 0.08289439062355086, 0.01715690310811624, -0.07257122172450181, -0.013804640759190079, 0.04770259846554836, 0.008318820257045445, -0.023639545008336427, -0.0037765424071949383, 0.00878780496532272, 0.001286487733523245, -0.00241815827996561, -0.00032699767544386305, 0.00047761004627489, 5.577576455095823e-05, -7.249498847983205e-05, -1.772932723653753e-05, -1.061598110485562e-05, -2.727912312749936e-05, -4.84062636949556e-05, -8.56389762576848e-05, -0.0001609044043018315, -0.00032300748651784006, -0.0007014603337645309, -0.0019198857946925232, -0.009595223127881816, 1.028978700851531]
        val = horner(coeffs, 400.000000000363798*(T - 647.092499999999973))
    elif 647.095 < T <= 647.0959:
        # (3.157060421406988e-13, 2.45822459571187e-13, 1.085522247937845e-12)
        coeffs = [-4.6759843826293945e-05, -5.772896111011505e-05, 0.00029446277767419815, 0.0003508694935590029, -0.0008547466713935137, -0.0009755055652931333, 0.0014907524455338717, 0.001618919734028168, -0.0017289135430473834, -0.0017798669869080186, 0.0013886251217627432, 0.0013546462596423225, -0.0007873586769164831, -0.0007338769426041836, 0.0003066651955805355, 0.0002761181634696186, -8.79955611026162e-05, -8.321968066127283e-05, 2.5915467887216437e-06, -1.7553803992598205e-06, -2.3491923775154078e-05, -3.44673285415098e-05, -4.9372020503035685e-05, -7.664722940048385e-05, -0.0001250530805668726, -0.00021831841691788423, -0.00042691414370452097, -0.001040188531514219, -0.004763499073651319, 1.0120628799127538]
        val = horner(coeffs, 2222.22222221862921*(T - 647.095450000000028))
    elif 647.0959 < T <= 647.09599:
        # 1e-12 max
        coeffs = [-0.001290641725063324, -0.0003778710961341858, 0.010535558685660362, 0.00301162526011467, -0.03902981849387288, -0.010885929688811302, 0.08678736770525575, 0.02358410635497421, -0.12912050032173283, -0.034109388157958165, 0.13565029118035454, 0.03472165907442104, -0.10352977068032487, -0.02556207562884083, 0.05817742309318419, 0.013771163636647543, -0.024138506512485947, -0.005434921958681116, 0.007345041980045153, 0.0015556455067553543, -0.0016143303365225847, -0.00032046819030284723, 0.0002439657241879445, 3.995263912059954e-05, -3.2193630974397536e-05, -1.4593553628650113e-05, -1.416486926064664e-05, -2.4128956016555847e-05, -3.973400246815029e-05, -6.928431227809506e-05, -0.0001353896286368922, -0.0003305923160938848, -0.0015796620258456111, 1.0039108886171622]
        val = horner(coeffs, 22222.2222221862903*(T - 647.095945000000029))
    elif 647.09599 < T <= 647.095999:
        # (3.667783433807127e-12, 2.9131002083651304e-12, 1.512957262282052e-11)
        coeffs = [-4.733548848889768e-06, -7.60948023525998e-06, 1.7879367078421637e-05, 3.080385795328766e-05, -3.0345471259352053e-05, -5.6281653996848036e-05, 2.7725049449145445e-05, 5.695053056342658e-05, -1.5959582896130087e-05, -3.672514257857529e-05, 3.8304293354940455e-06, 1.2481707166500655e-05, -3.0915939763787037e-06, -6.661048546519055e-06, -4.9108654751606196e-06, -7.2591069035787825e-06, -1.2581070591322074e-05, -2.1996436737614644e-05, -4.2936997922793374e-05, -0.00010488089397681088, -0.0005081514985303057, 1.0012486252556012]
        val = horner(coeffs, 222222.222783279401*(T - 647.095994499999961))
    elif 647.095999 < T <= 647.09599999:
#         (1.301907916559972e-11, 1.5329272639853585e-11, 9.685327631780702e-11)
        coeffs = [-342.1649169921875, -286.70733642578125, 3696.8092041015625, 3062.9874267578125, -18628.276138305664, -15248.527671813965, 58141.86198616028, 46974.429047584534, -125898.94741535187, -100287.08330655098, 200761.28558278084, 157482.6384627223, -244240.7832775116, -188416.60116776824, 231738.97573629767, 175548.21367172152, -173839.0422554016, -129094.15901541244, 103905.44825894013, 75496.35919268127, -49644.42389814614, -35215.162472276075, 18947.322668814828, 13087.622056136344, -5751.365423334828, -3856.5318880818013, 1377.1159714494315, 892.9758781142154, -256.84921118732746, -160.2506746801855, 36.651197303384606, 21.84684986680145, -3.9021237365971047, -2.1983816265822327, 0.2992404034464471, 0.15648404704313634, -0.01571472634215354, -0.007373118351861052, 0.0005188793214245813, 0.00019840994932129874, -1.665223554836448e-05, -1.1892982883288106e-05, -1.36274066943054e-05, -2.2384801104155527e-05, -4.5679905172298085e-05, -0.0001857746544898474, 1.0003799295366436]
        val = horner(coeffs, 2020202.04154185997*(T - 647.095999495000001))
    elif T <= Tc:
        # Gotta go to linear interp
        val = 1.0000546416597242 - 5.464165972424162e-05*(T-647.09599999)/(647.096-647.09599999)
    return val*rhoc

def iapws95_rho_err(rho, T, P_spec):
    tau = Tc/T
    delta = rho*rhoc_inv
    dAddelta_res_val = iapws95_dAr_ddelta(tau, delta)
    d2Ad2delta_res_val = iapws95_d2Ar_ddelta2(tau, delta)
    P_calc = (1.0 + dAddelta_res_val*delta)*rho*R95*T
    err = P_calc - P_spec
    derr = R95*T*(rho*(rho*d2Ad2delta_res_val + 644.0*dAddelta_res_val)
                + 103684.0)*9.644689633887581e-06 # 1/322**2
    return err, derr

def iapws97_rho_extrapolated(T, P):
    # Intended to extend the range using first derivatives
    # for use in iapws-95 solver.
    try:
        rho = iapws97_rho(T, P)
    except:
        if T > 2273.15 and P < 50E6:
            T_border = 2273.15
            Pref = 1e6
            Tref = 1e3
            P_inv = 1.0/P
            pi = P*1e-6
            tau = Tref/T_border
            # region 5 upwards T extrapolate drho_dT
            dGr_dpi = iapws97_dGr_dpi_region5(tau, pi) 
            dG_dpi = 1e6*P_inv + dGr_dpi
            rho = P/(R97*T_border*pi*dG_dpi)
            d2Gr_dpidtau = iapws97_d2Gr_dpidtau_region5(tau, pi) 
            x0 = (dGr_dpi + Pref*P_inv)
            x1 = R97*T_border*T_border
            drho_dT = -Pref/(x1*x0) + Pref*Tref*d2Gr_dpidtau/(x1*T_border*x0*x0)
            rho = rho + drho_dT*(T - T_border)
        elif T > 1073.15 and 50e6 <= P <= 100e6:
            T_border = 1073.15
            Pref = 1e6
            Tref = 540.0
            P_inv = 1.0/P
            pi = P*1e-6
            tau = Tref/T_border
            # region 5 upwards T extrapolate drho_dT
            dGr_dpi = iapws97_dGr_dpi_region2(tau, pi) 
            dG_dpi = 1e6*P_inv + dGr_dpi
            rho = P/(R97*T_border*pi*dG_dpi)
            d2Gr_dpidtau = iapws97_d2Gr_dpidtau_region2(tau, pi) 
            x0 = (dGr_dpi + Pref*P_inv)
            x1 = R97*T_border*T_border
            drho_dT = -Pref/(x1*x0) + Pref*Tref*d2Gr_dpidtau/(x1*T_border*x0*x0)
            drho = drho_dT*(T - T_border)
            if (rho + drho) > .1*rho:
                # Do not take the extrapolation if density has decreased too much
                rho = rho + drho
        elif T < 273.15:
            # region 1 - fairly incompressible, don't bother extrapolating
            if P > 100e6: P = 100e6
            rho = iapws97_region1_rho(273.15, P)
        else:
            # Don't bother extrapolating to higher P at this point, the density
            # change is small there
            if P > 100e6:
                P = 100e6
            if T > 1073.15:
                T = 1073.15
            rho = iapws97_rho(T, P)
    return rho

def iapws95_P(T, rho):
    tau = Tc/T
    delta = rho*rhoc_inv
    dAddelta_val = iapws95_dAr_ddelta(tau, delta) + 1.0/delta
    return (dAddelta_val*delta)*rho*R95*T

def iapws95_P_err(T, rho, P_spec):
    tau = Tc/T
    delta = rho*rhoc_inv
    dAddelta_val = iapws95_dAr_ddelta(tau, delta) + 1.0/delta
    err = (dAddelta_val*delta)*rho*R95*T - P_spec
    dP_dT = rho*R95*delta*(dAddelta_val - tau*iapws95_d2Ar_ddeltadtau(tau, delta))
    return err, dP_dT

def iapws_T(P, rho):
    try:
        T = iapws97_T(P, rho)
    except:
        pass


def iapws95_rho(T, P):
    MAX_RHO_STEP = 200.0 # iapws95_rho(250, 1e9) is a good point showing the advantage of this
    rho = iapws97_rho_extrapolated(T, P)
#    print(rho)
    # newton solver overhead is huge.
#            print(rho, 'initial guess', drho_dT, 'drho_dT')
            
    err, derr = iapws95_rho_err(rho, T, P)
    drho = - err/derr
    if drho < -MAX_RHO_STEP:
        drho = -MAX_RHO_STEP
    elif drho > MAX_RHO_STEP:
        drho = MAX_RHO_STEP
    rho_old = rho + drho
#    print(rho_old)
    
    err, derr = iapws95_rho_err(rho_old, T, P)
    drho = - err/derr
    if drho < -MAX_RHO_STEP:
        drho = -MAX_RHO_STEP
    elif drho > MAX_RHO_STEP:
        drho = MAX_RHO_STEP
    rho = rho_old + drho
    # Adding iterations check did not slow anything down.
    iterations = 2
    while (abs(rho_old - rho) > abs(1e-11*rho)) and iterations < 100:
        rho_old = rho
        err, derr = iapws95_rho_err(rho, T, P)
        drho = - err/derr
        if drho < -MAX_RHO_STEP:
            drho = -MAX_RHO_STEP
        elif drho > MAX_RHO_STEP:
            drho = MAX_RHO_STEP
        rho = rho + drho
        iterations += 1
#        print(rho, err)
    return rho

