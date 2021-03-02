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
"""

import numpy as np
import pytest
import chemicals
from chemicals.air import *
from math import *
from fluids.numerics import assert_close, assert_close1d, assert_close2d, linspace, logspace, derivative

from chemicals.air import TAU_MAX_EXP_87


def func_vs_naive_tester(func, func_naive, T_min=1.0, T_max=5000.0, rho_min=1e-5, rho_max=50000.0, N=400, Tc=132.6312, rhoc=10447.7):
    errs = []
    rerr = 0
    Ts = linspace(T_min,  T_max, N)
    rhoc_inv = 1.0/rhoc
    for i, T in enumerate(Ts):
        tau = Tc/T
        rhos = logspace(log10(rho_min), log10(rho_max), N)
        for rho in rhos:
            delta = rho*rhoc_inv
            val = func(tau, delta)
            val_naive = func_naive(tau, delta)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    AARD, std, max_err = rerr/N**2, np.std(errs), np.max(errs)
    return AARD, std, max_err




def test_lemmon2000_A0():
    assert_close(lemmon2000_air_A0(0.36842, 0.5), -17.026123512818458, rtol=1e-15)

    # if statement test
    assert_close(lemmon2000_air_A0(TAU_MAX_EXP_87*(1-1e-15), .5),
                 lemmon2000_air_A0(TAU_MAX_EXP_87*(1+1e-15), .5), rtol=1e-15)

    # Value computed with SymPy
    assert_close(lemmon2000_air_A0(300, .5), -3.54214522586700104, rtol=5e-15)


    # Check points that may fail
    lemmon2000_air_A0(1e20, 1e20)
    lemmon2000_air_A0(1e-10, 1e-10)



def test_lemmon2000_air_dA0_dtau():
    assert_close(derivative(lambda tau: lemmon2000_air_A0(tau, 13000/11183.9), 126.192/200.0, dx=1e-7),
                 lemmon2000_air_dA0_dtau(126.192/200.0, 13000/11183.9), rtol=1e-9)

    assert_close(lemmon2000_air_dA0_dtau(0.36842, 0.5), 6.764336610288353, rtol=1e-14)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_A0(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_dA0_dtau(rat, .5))


def test_lemmon2000_air_d2A0_dtau2():
    assert_close(lemmon2000_air_d2A0_dtau2(126.192/200.0, 13000/11183.9), -6.260482844274295, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_dA0_dtau(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d2A0_dtau2(rat, .5))


def test_lemmon2000_air_d3A0_dtau3():
    assert_close(lemmon2000_air_d3A0_dtau3(0.36842, .5), 102.9144884392338, rtol=2e-12)
    assert_close(derivative(lambda tau: lemmon2000_air_d2A0_dtau2(tau, .5), 0.36842, dx=1e-7),
                 lemmon2000_air_d3A0_dtau3(0.36842, .5))

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d2A0_dtau2(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d3A0_dtau3(rat, .5))

def test_lemmon2000_air_d4A0_dtau4():
    assert_close(lemmon2000_air_d4A0_dtau4(126.192/200.0, 13000/11183.9), -94.8155327278803, rtol=1e-13)

    assert_close(derivative(lambda tau: lemmon2000_air_d3A0_dtau3(tau, .5), 0.36842, dx=4e-7),
                 lemmon2000_air_d4A0_dtau4(0.36842, .5))

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0,2.5, 2.0,  1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d3A0_dtau3(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d4A0_dtau4(rat, .5))


# Naive functions
def lemmon2000_air_Ar_naive(tau, delta):
    return 0.000233594806141999996*delta**11*tau**3.25*exp(-delta**2) - 0.0122523554252999996*delta**6*tau**1.25*exp(-delta) + 0.000164957183186000006*delta**6*tau**1.35000000000000009 - 0.0472103183731000034*delta**5*tau**0.949999999999999956*exp(-delta) - 0.042053322884200002*delta**4*tau**0.200000000000000011 + 0.0349008431981999989*delta**4*tau**0.349999999999999978 + 0.0112626704218000001*delta**4 + 0.134211176704000013*delta**3*tau**0.149999999999999994 - 0.17381369096999999*delta**3*tau**0.800000000000000044*exp(-delta) - 0.00938782884667000057*delta**3*tau**15*exp(-delta**3) - 0.031605587982100003*delta**3*tau**6*exp(-delta**2) - 0.0865421396646000041*delta**3 + 0.0714140178971000017*delta**2 + 0.713116392079000017*delta*tau**0.330000000000000016 - 1.61824192067000006*delta*tau**1.01000000000000001 - 0.101365037911999994*delta*tau**1.60000000000000009*exp(-delta) + 0.0148287891978000005*delta*tau**3.5*exp(-delta**3) - 0.146629609712999986*delta*tau**3.60000000000000009*exp(-delta**2) + 0.118160747228999996*delta


def lemmon2000_air_dAr_dtau_naive(tau, delta):
    return 0.000759183119961500015*delta**11*tau**2.25*exp(-delta**2) - 0.0153154442816249986*delta**6*tau**0.25*exp(-delta) + 0.00022269219730110002*delta**6*tau**0.350000000000000089 - 0.0448498024544450036*delta**5*tau**(-0.0500000000000000444)*exp(-delta) - 0.0084106645768400011*delta**4*tau**(-0.800000000000000044) + 0.0122152951193699993*delta**4*tau**(-0.650000000000000022) + 0.0201316765056000005*delta**3*tau**(-0.849999999999999978) - 0.139050952775999992*delta**3*tau**(-0.199999999999999956)*exp(-delta) - 0.14081743270005001*delta**3*tau**14*exp(-delta**3) - 0.189633527892600018*delta**3*tau**5*exp(-delta**2) + 0.235328409386070025*delta*tau**(-0.669999999999999929) - 1.63442433987669999*delta*tau**0.0100000000000000089 - 0.162184060659199991*delta*tau**0.600000000000000089*exp(-delta) + 0.0519007621922999984*delta*tau**2.5*exp(-delta**3) - 0.527866594966799996*delta*tau**2.60000000000000009*exp(-delta**2)

def lemmon2000_air_d2Ar_dtau2_naive(tau, delta):
    return delta*(0.00170816201991337503*delta**10*tau**1.25*exp(-delta**2) - 0.00382886107040624965*delta**5*tau**(-0.75)*exp(-delta) + 0.0000779422690553850206*delta**5*tau**(-0.649999999999999911) + 0.00224249012272225226*delta**4*tau**(-1.05000000000000004)*exp(-delta) + 0.00672853166147200157*delta**3*tau**(-1.80000000000000004) - 0.00793994182759050031*delta**3*tau**(-1.64999999999999991) - 0.0171119250297600001*delta**2*tau**(-1.85000000000000009) + 0.0278101905551999921*delta**2*tau**(-1.19999999999999996)*exp(-delta) - 1.9714440578007002*delta**2*tau**13*exp(-delta**3) - 0.948167639463000089*delta**2*tau**4*exp(-delta**2) - 0.15767003428866691*tau**(-1.66999999999999993) - 0.0163442433987670138*tau**(-0.989999999999999991) - 0.0973104363955200058*tau**(-0.399999999999999911)*exp(-delta) + 0.129751905480749996*tau**1.5*exp(-delta**3) - 1.37245314691368003*tau**1.60000000000000009*exp(-delta**2))

def lemmon2000_air_d3Ar_dtau3_naive(tau, delta):
    return delta*(0.00213520252489171874*delta**10*tau**0.25*exp(-delta**2) + 0.00287164580280468724*delta**5*tau**(-1.75)*exp(-delta) - 0.000050662474886000258*delta**5*tau**(-1.64999999999999991) - 0.00235461462885836479*delta**4*tau**(-2.04999999999999982)*exp(-delta) - 0.0121113569906496025*delta**3*tau**(-2.79999999999999982) + 0.0131009040155243249*delta**3*tau**(-2.64999999999999991) + 0.0316570613050560015*delta**2*tau**(-2.85000000000000009) - 0.0333722286662399906*delta**2*tau**(-2.20000000000000018)*exp(-delta) - 25.6287727514091017*delta**2*tau**12*exp(-delta**3) - 3.79267055785200036*delta**2*tau**3*exp(-delta**2) + 0.263308957262073706*tau**(-2.66999999999999993) + 0.0161808009647793419*tau**(-1.98999999999999999) + 0.0389241745582079926*tau**(-1.39999999999999991)*exp(-delta) + 0.194627858221124994*tau**0.5*exp(-delta**3) - 2.19592503506188796*tau**0.600000000000000089*exp(-delta**2))

def lemmon2000_air_d4Ar_dtau4_naive(tau, delta):
    return delta*(0.000533800631222929685*delta**10*tau**(-0.75)*exp(-delta**2) - 0.00502538015490820288*delta**5*tau**(-2.75)*exp(-delta) + 0.0000835930835619004249*delta**5*tau**(-2.64999999999999991) + 0.00482695998915964701*delta**4*tau**(-3.04999999999999982)*exp(-delta) + 0.0339117995738188877*delta**3*tau**(-3.79999999999999982) - 0.0347173956411394591*delta**3*tau**(-3.64999999999999991) - 0.0902226247194096026*delta**2*tau**(-3.85000000000000009) + 0.0734189030657279862*delta**2*tau**(-3.20000000000000018)*exp(-delta) - 307.545273016909221*delta**2*tau**11*exp(-delta**3) - 11.3780116735560011*delta**2*tau**2*exp(-delta**2) - 0.703034915889736767*tau**(-3.66999999999999993) - 0.0321997939199108879*tau**(-2.99000000000000021) - 0.0544938443814911855*tau**(-2.39999999999999991)*exp(-delta) + 0.0973139291105624971*tau**(-0.5)*exp(-delta**3) - 1.31755502103713296*tau**(-0.399999999999999911)*exp(-delta**2))


def lemmon2000_air_dAr_ddelta_naive(tau, delta):
    return -0.000467189612283999993*delta**12*tau**3.25*exp(-delta**2) + 0.00256954286756199985*delta**10*tau**3.25*exp(-delta**2) + 0.0122523554252999996*delta**6*tau**1.25*exp(-delta) + 0.0472103183731000034*delta**5*tau**0.949999999999999956*exp(-delta) + 0.0281634865400100035*delta**5*tau**15*exp(-delta**3) - 0.0735141325518000044*delta**5*tau**1.25*exp(-delta) + 0.000989743099116000089*delta**5*tau**1.35000000000000009 - 0.23605159186550001*delta**4*tau**0.949999999999999956*exp(-delta) + 0.063211175964200006*delta**4*tau**6*exp(-delta**2) - 0.168213291536800008*delta**3*tau**0.200000000000000011 + 0.139603372792799996*delta**3*tau**0.349999999999999978 + 0.17381369096999999*delta**3*tau**0.800000000000000044*exp(-delta) - 0.0444863675934000016*delta**3*tau**3.5*exp(-delta**3) + 0.0450506816872000004*delta**3 + 0.402633530112000038*delta**2*tau**0.149999999999999994 - 0.52144107290999997*delta**2*tau**0.800000000000000044*exp(-delta) - 0.0281634865400100035*delta**2*tau**15*exp(-delta**3) - 0.0948167639463000089*delta**2*tau**6*exp(-delta**2) + 0.293259219425999973*delta**2*tau**3.60000000000000009*exp(-delta**2) - 0.259626418993800012*delta**2 + 0.101365037911999994*delta*tau**1.60000000000000009*exp(-delta) + 0.142828035794200003*delta + 0.713116392079000017*tau**0.330000000000000016 - 1.61824192067000006*tau**1.01000000000000001 - 0.101365037911999994*tau**1.60000000000000009*exp(-delta) + 0.0148287891978000005*tau**3.5*exp(-delta**3) - 0.146629609712999986*tau**3.60000000000000009*exp(-delta**2) + 0.118160747228999996

def lemmon2000_air_d2Ar_ddelta2_naive(tau, delta):
    return 0.000934379224567999985*delta**13*tau**3.25*exp(-delta**2) - 0.0107453610825320005*delta**11*tau**3.25*exp(-delta**2) + 0.0256954286756199968*delta**9*tau**3.25*exp(-delta**2) - 0.0844904596200300173*delta**7*tau**15*exp(-delta**3) - 0.0122523554252999996*delta**6*tau**1.25*exp(-delta) - 0.0472103183731000034*delta**5*tau**0.949999999999999956*exp(-delta) - 0.126422351928400012*delta**5*tau**6*exp(-delta**2) + 0.147028265103600009*delta**5*tau**1.25*exp(-delta) + 0.133459102780200012*delta**5*tau**3.5*exp(-delta**3) + 0.47210318373100002*delta**4*tau**0.949999999999999956*exp(-delta) + 0.225307892320080028*delta**4*tau**15*exp(-delta**3) - 0.367570662759000022*delta**4*tau**1.25*exp(-delta) + 0.00494871549558000001*delta**4*tau**1.35000000000000009 - 0.17381369096999999*delta**3*tau**0.800000000000000044*exp(-delta) - 0.94420636746200004*delta**3*tau**0.949999999999999956*exp(-delta) + 0.442478231749400042*delta**3*tau**6*exp(-delta**2) - 0.586518438851999946*delta**3*tau**3.60000000000000009*exp(-delta**2) - 0.504639874610400052*delta**2*tau**0.200000000000000011 + 0.418810118378399987*delta**2*tau**0.349999999999999978 + 1.04288214581999994*delta**2*tau**0.800000000000000044*exp(-delta) - 0.177945470373600007*delta**2*tau**3.5*exp(-delta**3) + 0.135152045061600001*delta**2 + 0.805267060224000075*delta*tau**0.149999999999999994 - 1.04288214581999994*delta*tau**0.800000000000000044*exp(-delta) - 0.0563269730800200069*delta*tau**15*exp(-delta**3) - 0.189633527892600018*delta*tau**6*exp(-delta**2) - 0.101365037911999994*delta*tau**1.60000000000000009*exp(-delta) + 0.879777658277999919*delta*tau**3.60000000000000009*exp(-delta**2) - 0.519252837987600024*delta + 0.202730075823999989*tau**1.60000000000000009*exp(-delta) + 0.142828035794200003

def lemmon2000_air_d3Ar_ddelta3_naive(tau, delta):
    return -0.00186875844913599997*delta**14*tau**3.25*exp(-delta**2) + 0.0336376520844480012*delta**12*tau**3.25*exp(-delta**2) - 0.169589829259091995*delta**10*tau**3.25*exp(-delta**2) + 0.25347137886009008*delta**9*tau**15*exp(-delta**3) + 0.231258858080579971*delta**8*tau**3.25*exp(-delta**2) - 0.400377308340600035*delta**7*tau**3.5*exp(-delta**3) - 1.26735689430045029*delta**6*tau**15*exp(-delta**3) + 0.252844703856800024*delta**6*tau**6*exp(-delta**2) + 0.0122523554252999996*delta**6*tau**1.25*exp(-delta) + 0.0472103183731000034*delta**5*tau**0.949999999999999956*exp(-delta) - 0.220542397655400013*delta**5*tau**1.25*exp(-delta) - 0.708154775596500086*delta**4*tau**0.949999999999999956*exp(-delta) - 1.51706822314080014*delta**4*tau**6*exp(-delta**2) + 1.10271198827700001*delta**4*tau**1.25*exp(-delta) + 1.20113192502180022*delta**4*tau**3.5*exp(-delta**3) + 1.17303687770399989*delta**4*tau**3.60000000000000009*exp(-delta**2) + 0.17381369096999999*delta**3*tau**0.800000000000000044*exp(-delta) + 2.83261910238600034*delta**3*tau**0.949999999999999956*exp(-delta) + 1.07021248852038009*delta**3*tau**15*exp(-delta**3) - 1.47028265103600009*delta**3*tau**1.25*exp(-delta) + 0.01979486198232*delta**3*tau**1.35000000000000009 - 1.5643232187299998*delta**2*tau**0.800000000000000044*exp(-delta) - 2.83261910238600034*delta**2*tau**0.949999999999999956*exp(-delta) + 1.70670175103340016*delta**2*tau**6*exp(-delta**2) - 3.51911063311199968*delta**2*tau**3.60000000000000009*exp(-delta**2) - 1.0092797492208001*delta*tau**0.200000000000000011 + 0.837620236756799974*delta*tau**0.349999999999999978 + 3.1286464374599996*delta*tau**0.800000000000000044*exp(-delta) + 0.101365037911999994*delta*tau**1.60000000000000009*exp(-delta) - 0.355890940747200013*delta*tau**3.5*exp(-delta**3) + 0.270304090123200003*delta + 0.805267060224000075*tau**0.149999999999999994 - 1.04288214581999994*tau**0.800000000000000044*exp(-delta) - 0.0563269730800200069*tau**15*exp(-delta**3) - 0.189633527892600018*tau**6*exp(-delta**2) - 0.304095113735999956*tau**1.60000000000000009*exp(-delta) + 0.879777658277999919*tau**3.60000000000000009*exp(-delta**2) - 0.519252837987600024

def lemmon2000_air_d4Ar_ddelta4_naive(tau, delta):
    return 0.00373751689827199994*delta**15*tau**3.25*exp(-delta**2) - 0.0934379224567999933*delta**13*tau**3.25*exp(-delta**2) - 0.760414136580270239*delta**11*tau**15*exp(-delta**3) + 0.742831483531560033*delta**11*tau**3.25*exp(-delta**2) - 2.15841600875207984*delta**9*tau**3.25*exp(-delta**2) + 1.20113192502180022*delta**9*tau**3.5*exp(-delta**3) + 6.08331309264216191*delta**8*tau**15*exp(-delta**3) - 0.505689407713600048*delta**7*tau**6*exp(-delta**2) + 1.85007086464463977*delta**7*tau**3.25*exp(-delta**2) - 0.0122523554252999996*delta**6*tau**1.25*exp(-delta) - 6.40603693344960057*delta**6*tau**3.5*exp(-delta**3) - 0.0472103183731000034*delta**5*tau**0.949999999999999956*exp(-delta) - 10.8147788313638422*delta**5*tau**15*exp(-delta**3) + 4.55120466942240043*delta**5*tau**6*exp(-delta**2) + 0.294056530207200018*delta**5*tau**1.25*exp(-delta) - 2.34607375540799978*delta**5*tau**3.60000000000000009*exp(-delta**2) + 0.94420636746200004*delta**4*tau**0.949999999999999956*exp(-delta) - 2.20542397655400002*delta**4*tau**1.25*exp(-delta) - 0.17381369096999999*delta**3*tau**0.800000000000000044*exp(-delta) - 5.66523820477200069*delta**3*tau**0.949999999999999956*exp(-delta) - 9.48167639463*delta**3*tau**6*exp(-delta**2) + 5.88113060414400035*delta**3*tau**1.25*exp(-delta) + 5.87220052232880096*delta**3*tau**3.5*exp(-delta**3) + 11.7303687770399989*delta**3*tau**3.60000000000000009*exp(-delta**2) + 2.08576429163999988*delta**2*tau**0.800000000000000044*exp(-delta) + 11.3304764095440014*delta**2*tau**0.949999999999999956*exp(-delta) + 3.37961838480120003*delta**2*tau**15*exp(-delta**3) - 4.41084795310800004*delta**2*tau**1.25*exp(-delta) + 0.0593845859469600001*delta**2*tau**1.35000000000000009 - 6.25729287491999919*delta*tau**0.800000000000000044*exp(-delta) - 5.66523820477200069*delta*tau**0.949999999999999956*exp(-delta) + 3.79267055785200036*delta*tau**6*exp(-delta**2) - 0.101365037911999994*delta*tau**1.60000000000000009*exp(-delta) - 8.79777658277999919*delta*tau**3.60000000000000009*exp(-delta**2) - 1.0092797492208001*tau**0.200000000000000011 + 0.837620236756799974*tau**0.349999999999999978 + 4.17152858327999976*tau**0.800000000000000044*exp(-delta) + 0.405460151647999978*tau**1.60000000000000009*exp(-delta) - 0.355890940747200013*tau**3.5*exp(-delta**3) + 0.270304090123200003


def lemmon2000_air_d2Ar_ddeltadtau_naive(tau, delta):
    return -0.00151836623992300003*delta**12*tau**2.25*exp(-delta**2) + 0.0083510143195764993*delta**10*tau**2.25*exp(-delta**2) + 0.0153154442816249986*delta**6*tau**0.25*exp(-delta) + 0.0448498024544450036*delta**5*tau**(-0.0500000000000000444)*exp(-delta) - 0.0918926656897500055*delta**5*tau**0.25*exp(-delta) + 0.00133615318380660023*delta**5*tau**0.350000000000000089 + 0.422452298100150059*delta**5*tau**14*exp(-delta**3) - 0.22424901227222499*delta**4*tau**(-0.0500000000000000444)*exp(-delta) + 0.379267055785200036*delta**4*tau**5*exp(-delta**2) - 0.0336426583073600044*delta**3*tau**(-0.800000000000000044) + 0.0488611804774799971*delta**3*tau**(-0.650000000000000022) + 0.139050952775999992*delta**3*tau**(-0.199999999999999956)*exp(-delta) - 0.155702286576899995*delta**3*tau**2.5*exp(-delta**3) + 0.0603950295168000015*delta**2*tau**(-0.849999999999999978) - 0.417152858327999976*delta**2*tau**(-0.199999999999999956)*exp(-delta) - 0.422452298100150059*delta**2*tau**14*exp(-delta**3) - 0.568900583677800054*delta**2*tau**5*exp(-delta**2) + 1.05573318993359999*delta**2*tau**2.60000000000000009*exp(-delta**2) + 0.162184060659199991*delta*tau**0.600000000000000089*exp(-delta) + 0.235328409386070025*tau**(-0.669999999999999929) - 1.63442433987669999*tau**0.0100000000000000089 - 0.162184060659199991*tau**0.600000000000000089*exp(-delta) + 0.0519007621922999984*tau**2.5*exp(-delta**3) - 0.527866594966799996*tau**2.60000000000000009*exp(-delta**2)

def lemmon2000_air_d3Ar_ddeltadtau2_naive(tau, delta):
    return -0.00341632403982675007*delta**12*tau**1.25*exp(-delta**2) + 0.0187897822190471221*delta**10*tau**1.25*exp(-delta**2) + 0.00382886107040624965*delta**6*tau**(-0.75)*exp(-delta) - 0.00224249012272225226*delta**5*tau**(-1.05000000000000004)*exp(-delta) - 0.0229731664224375014*delta**5*tau**(-0.75)*exp(-delta) + 0.000467653614332310178*delta**5*tau**(-0.649999999999999911) + 5.9143321734021006*delta**5*tau**13*exp(-delta**3) + 0.0112124506136112596*delta**4*tau**(-1.05000000000000004)*exp(-delta) + 1.89633527892600018*delta**4*tau**4*exp(-delta**2) + 0.0269141266458880063*delta**3*tau**(-1.80000000000000004) - 0.0317597673103620012*delta**3*tau**(-1.64999999999999991) - 0.0278101905551999921*delta**3*tau**(-1.19999999999999996)*exp(-delta) - 0.389255716442249988*delta**3*tau**1.5*exp(-delta**3) - 0.0513357750892800002*delta**2*tau**(-1.85000000000000009) + 0.083430571665599973*delta**2*tau**(-1.19999999999999996)*exp(-delta) - 5.9143321734021006*delta**2*tau**13*exp(-delta**3) - 2.84450291838900027*delta**2*tau**4*exp(-delta**2) + 2.74490629382736007*delta**2*tau**1.60000000000000009*exp(-delta**2) + 0.0973104363955200058*delta*tau**(-0.399999999999999911)*exp(-delta) - 0.15767003428866691*tau**(-1.66999999999999993) - 0.0163442433987670138*tau**(-0.989999999999999991) - 0.0973104363955200058*tau**(-0.399999999999999911)*exp(-delta) + 0.129751905480749996*tau**1.5*exp(-delta**3) - 1.37245314691368003*tau**1.60000000000000009*exp(-delta**2)

def lemmon2000_air_d3Ar_ddelta2dtau_naive(tau, delta):
    return 0.00303673247984600006*delta**13*tau**2.25*exp(-delta**2) - 0.0349224235182290024*delta**11*tau**2.25*exp(-delta**2) + 0.0835101431957649964*delta**9*tau**2.25*exp(-delta**2) - 1.26735689430045029*delta**7*tau**14*exp(-delta**3) - 0.0153154442816249986*delta**6*tau**0.25*exp(-delta) - 0.0448498024544450036*delta**5*tau**(-0.0500000000000000444)*exp(-delta) + 0.183785331379500011*delta**5*tau**0.25*exp(-delta) - 0.758534111570400071*delta**5*tau**5*exp(-delta**2) + 0.467106859730700041*delta**5*tau**2.5*exp(-delta**3) + 0.44849802454444998*delta**4*tau**(-0.0500000000000000444)*exp(-delta) - 0.459463328448750041*delta**4*tau**0.25*exp(-delta) + 0.00668076591903300071*delta**4*tau**0.350000000000000089 + 3.37961838480120047*delta**4*tau**14*exp(-delta**3) - 0.139050952775999992*delta**3*tau**(-0.199999999999999956)*exp(-delta) - 0.89699604908889996*delta**3*tau**(-0.0500000000000000444)*exp(-delta) + 2.65486939049640025*delta**3*tau**5*exp(-delta**2) - 2.11146637986719998*delta**3*tau**2.60000000000000009*exp(-delta**2) - 0.100927974922080013*delta**2*tau**(-0.800000000000000044) + 0.146583541432439984*delta**2*tau**(-0.650000000000000022) + 0.834305716655999952*delta**2*tau**(-0.199999999999999956)*exp(-delta) - 0.622809146307599981*delta**2*tau**2.5*exp(-delta**3) + 0.120790059033600003*delta*tau**(-0.849999999999999978) - 0.834305716655999952*delta*tau**(-0.199999999999999956)*exp(-delta) - 0.162184060659199991*delta*tau**0.600000000000000089*exp(-delta) - 0.844904596200300118*delta*tau**14*exp(-delta**3) - 1.13780116735560011*delta*tau**5*exp(-delta**2) + 3.16719956980079997*delta*tau**2.60000000000000009*exp(-delta**2) + 0.324368121318399982*tau**0.600000000000000089*exp(-delta)


def lemmon2000_air_d4Ar_ddelta2dtau2_naive(tau, delta):
    return 0.00683264807965350014*delta**13*tau**1.25*exp(-delta**2) - 0.0785754529160152537*delta**11*tau**1.25*exp(-delta**2) + 0.187897822190471242*delta**9*tau**1.25*exp(-delta**2) - 17.7429965202063045*delta**7*tau**13*exp(-delta**3) - 0.00382886107040624965*delta**6*tau**(-0.75)*exp(-delta) + 0.00224249012272225226*delta**5*tau**(-1.05000000000000004)*exp(-delta) + 0.0459463328448750027*delta**5*tau**(-0.75)*exp(-delta) - 3.79267055785200036*delta**5*tau**4*exp(-delta**2) + 1.16776714932675008*delta**5*tau**1.5*exp(-delta**3) - 0.0224249012272225191*delta**4*tau**(-1.05000000000000004)*exp(-delta) - 0.11486583211218751*delta**4*tau**(-0.75)*exp(-delta) + 0.00233826807166155094*delta**4*tau**(-0.649999999999999911) + 47.3146573872168048*delta**4*tau**13*exp(-delta**3) + 0.0278101905551999921*delta**3*tau**(-1.19999999999999996)*exp(-delta) + 0.0448498024544450383*delta**3*tau**(-1.05000000000000004)*exp(-delta) + 13.2743469524820021*delta**3*tau**4*exp(-delta**2) - 5.48981258765472013*delta**3*tau**1.60000000000000009*exp(-delta**2) + 0.0807423799376640189*delta**2*tau**(-1.80000000000000004) - 0.0952793019310859968*delta**2*tau**(-1.64999999999999991) - 0.166861143331199946*delta**2*tau**(-1.19999999999999996)*exp(-delta) - 1.55702286576899995*delta**2*tau**1.5*exp(-delta**3) - 0.10267155017856*delta*tau**(-1.85000000000000009) + 0.166861143331199946*delta*tau**(-1.19999999999999996)*exp(-delta) - 0.0973104363955200058*delta*tau**(-0.399999999999999911)*exp(-delta) - 11.8286643468042012*delta*tau**13*exp(-delta**3) - 5.68900583677800054*delta*tau**4*exp(-delta**2) + 8.23471888148208109*delta*tau**1.60000000000000009*exp(-delta**2) + 0.194620872791040012*tau**(-0.399999999999999911)*exp(-delta)

def lemmon2000_air_d4Ar_ddeltadtau3_naive(tau, delta):
    return -0.00427040504978343748*delta**12*tau**0.25*exp(-delta**2) + 0.0234872277738089018*delta**10*tau**0.25*exp(-delta**2) - 0.00287164580280468724*delta**6*tau**(-1.75)*exp(-delta) + 0.00235461462885836479*delta**5*tau**(-2.04999999999999982)*exp(-delta) + 0.0172298748168281252*delta**5*tau**(-1.75)*exp(-delta) - 0.000303974849316001575*delta**5*tau**(-1.64999999999999991) + 76.8863182542273051*delta**5*tau**12*exp(-delta**3) - 0.0117730731442918235*delta**4*tau**(-2.04999999999999982)*exp(-delta) + 7.58534111570400071*delta**4*tau**3*exp(-delta**2) - 0.0484454279625984099*delta**3*tau**(-2.79999999999999982) + 0.0524036160620972996*delta**3*tau**(-2.64999999999999991) + 0.0333722286662399906*delta**3*tau**(-2.20000000000000018)*exp(-delta) - 0.583883574663375038*delta**3*tau**0.5*exp(-delta**3) + 0.0949711839151680115*delta**2*tau**(-2.85000000000000009) - 0.100116685998719965*delta**2*tau**(-2.20000000000000018)*exp(-delta) + 4.39185007012377593*delta**2*tau**0.600000000000000089*exp(-delta**2) - 76.8863182542273051*delta**2*tau**12*exp(-delta**3) - 11.3780116735560011*delta**2*tau**3*exp(-delta**2) - 0.0389241745582079926*delta*tau**(-1.39999999999999991)*exp(-delta) + 0.263308957262073706*tau**(-2.66999999999999993) + 0.0161808009647793419*tau**(-1.98999999999999999) + 0.0389241745582079926*tau**(-1.39999999999999991)*exp(-delta) + 0.194627858221124994*tau**0.5*exp(-delta**3) - 2.19592503506188796*tau**0.600000000000000089*exp(-delta**2)

def lemmon2000_air_d4Ar_ddelta3dtau_naive(tau, delta):
    return -0.00607346495969200012*delta**14*tau**2.25*exp(-delta**2) + 0.109322369274456002*delta**12*tau**2.25*exp(-delta**2) - 0.551166945092048999*delta**10*tau**2.25*exp(-delta**2) + 3.80207068290135108*delta**9*tau**14*exp(-delta**3) + 0.751591288761884857*delta**8*tau**2.25*exp(-delta**2) - 1.40132057919210018*delta**7*tau**2.5*exp(-delta**3) + 0.0153154442816249986*delta**6*tau**0.25*exp(-delta) - 19.0103534145067528*delta**6*tau**14*exp(-delta**3) + 1.51706822314080014*delta**6*tau**5*exp(-delta**2) + 0.0448498024544450036*delta**5*tau**(-0.0500000000000000444)*exp(-delta) - 0.275677997069250003*delta**5*tau**0.25*exp(-delta) - 0.672747036816675026*delta**4*tau**(-0.0500000000000000444)*exp(-delta) + 1.37838998534625001*delta**4*tau**0.25*exp(-delta) - 9.10240933884480086*delta**4*tau**5*exp(-delta**2) + 4.20396173757630098*delta**4*tau**2.5*exp(-delta**3) + 4.22293275973439997*delta**4*tau**2.60000000000000009*exp(-delta**2) + 0.139050952775999992*delta**3*tau**(-0.199999999999999956)*exp(-delta) + 2.6909881472667001*delta**3*tau**(-0.0500000000000000444)*exp(-delta) - 1.83785331379500017*delta**3*tau**0.25*exp(-delta) + 0.0267230636761320028*delta**3*tau**0.350000000000000089 + 16.053187327805702*delta**3*tau**14*exp(-delta**3) - 1.25145857498399993*delta**2*tau**(-0.199999999999999956)*exp(-delta) - 2.6909881472667001*delta**2*tau**(-0.0500000000000000444)*exp(-delta) + 10.2402105062004019*delta**2*tau**5*exp(-delta**2) - 12.6687982792031999*delta**2*tau**2.60000000000000009*exp(-delta**2) - 0.201855949844160026*delta*tau**(-0.800000000000000044) + 0.293167082864879969*delta*tau**(-0.650000000000000022) + 2.50291714996799985*delta*tau**(-0.199999999999999956)*exp(-delta) + 0.162184060659199991*delta*tau**0.600000000000000089*exp(-delta) - 1.24561829261519996*delta*tau**2.5*exp(-delta**3) + 0.120790059033600003*tau**(-0.849999999999999978) - 0.834305716655999952*tau**(-0.199999999999999956)*exp(-delta) - 0.486552181977599973*tau**0.600000000000000089*exp(-delta) - 0.844904596200300118*tau**14*exp(-delta**3) - 1.13780116735560011*tau**5*exp(-delta**2) + 3.16719956980079997*tau**2.60000000000000009*exp(-delta**2)

def test_lemmon2000_air_Ar():
    assert_close(lemmon2000_air_Ar(0.36842, 0.15880050154579475), 0.004798812280624336, rtol=1e-13)

def test_lemmon2000_air_dAr_dtau():
    assert_close(lemmon2000_air_dAr_dtau(0.36842, 0.15880050154579475),  -0.20189573196786642, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_Ar(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_dAr_dtau(rat, .5))

def test_lemmon2000_air_d2Ar_dtau2():
    assert_close(lemmon2000_air_d2Ar_dtau2(132.6312/200.0, 13000/10447.7), -0.7632109061747537, rtol=1e-14)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_dAr_dtau(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d2Ar_dtau2(rat, .5))

def test_lemmon2000_air_d3Ar_dtau3():
    assert_close(lemmon2000_air_d3Ar_dtau3(132.6312/200.0, 13000/10447.7), 0.27922007457420017, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d2Ar_dtau2(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d3Ar_dtau3(rat, .5))

def test_lemmon2000_air_d4Ar_dtau4():
    assert_close(lemmon2000_air_d4Ar_dtau4(132.6312/200.0, 13000/10447.7), -8.197368061417675, rtol=1e-14)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d3Ar_dtau3(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d4Ar_dtau4(rat, .5))


def test_lemmon2000_air_dAr_ddelta():
    assert_close(lemmon2000_air_dAr_ddelta(0.36842, 0.15880050154579475), 0.0428706712678839, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda delta: lemmon2000_air_Ar(.5, delta), rat, dx=rat*1e-5),
                     lemmon2000_air_dAr_ddelta(.5, rat))

def test_lemmon2000_air_d2Ar_ddelta2():
    assert_close(lemmon2000_air_d2Ar_ddelta2(0.36842, 0.15880050154579475),  0.15165011962677075, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda delta: lemmon2000_air_dAr_ddelta(.5, delta), rat, dx=rat*1e-5),
                     lemmon2000_air_d2Ar_ddelta2(.5, rat))

def test_lemmon2000_air_d3Ar_ddelta3():
    assert_close(lemmon2000_air_d3Ar_ddelta3(0.36842, 0.15880050154579475), -0.09682239073996685, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda delta: lemmon2000_air_d2Ar_ddelta2(.5, delta), rat, dx=rat*1e-5),
                     lemmon2000_air_d3Ar_ddelta3(.5, rat))

def test_lemmon2000_air_d4Ar_ddelta4():
    assert_close(lemmon2000_air_d4Ar_ddelta4(0.36842, 0.15880050154579475),  1.0647113576834495, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda delta: lemmon2000_air_d3Ar_ddelta3(.5, delta), rat, dx=rat*1e-5),
                     lemmon2000_air_d4Ar_ddelta4(.5, rat))

def test_lemmon2000_air_d2Ar_ddeltadtau():
    assert_close(lemmon2000_air_d2Ar_ddeltadtau(0.36842, 0.15880050154579475), -1.261887383081615, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_dAr_ddelta(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d2Ar_ddeltadtau(rat, .5))

def test_lemmon2000_air_d3Ar_ddeltadtau2():
    assert_close(lemmon2000_air_d3Ar_ddeltadtau2(132.6312/200.0, 13000/10447.7), -0.19089212184849963, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d2Ar_ddeltadtau(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d3Ar_ddeltadtau2(rat, .5))

def test_lemmon2000_air_d3Ar_ddelta2dtau():
    assert_close(lemmon2000_air_d3Ar_ddelta2dtau(132.6312/200.0, 13000/10447.7), 0.014417881989408202, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d2Ar_ddelta2(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d3Ar_ddelta2dtau(rat, .5))

def test_lemmon2000_air_d3Ar_ddelta2dtau():
    assert_close(lemmon2000_air_d4Ar_ddelta2dtau2(132.6312/200.0, 13000/10447.7),  0.11968731127306076, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d3Ar_ddelta2dtau(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d4Ar_ddelta2dtau2(rat, .5))

def test_lemmon2000_air_d4Ar_ddeltadtau3():
    assert_close(lemmon2000_air_d4Ar_ddeltadtau3(132.6312/200.0, 13000/10447.7), 2.077739387492131, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d3Ar_ddeltadtau2(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d4Ar_ddeltadtau3(rat, .5))

def test_lemmon2000_air_d4Ar_ddelta3dtau():

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d3Ar_ddelta3(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d4Ar_ddelta3dtau(rat, .5))

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_Ar_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_Ar, lemmon2000_air_Ar_naive, N=100)
    assert AARD < 1e-13
    # If enough points happen, can find some pretty big discrepancie
    assert max_err < 1e-8
#    print(AARD, std, max_err)

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_dAr_dtau_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_dAr_dtau, lemmon2000_air_dAr_dtau_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8


@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d2Ar_dtau2_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d2Ar_dtau2, lemmon2000_air_d2Ar_dtau2_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d3Ar_dtau3_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d3Ar_dtau3, lemmon2000_air_d3Ar_dtau3_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d4Ar_dtau4_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d4Ar_dtau4, lemmon2000_air_d4Ar_dtau4_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_dAr_ddelta_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_dAr_ddelta, lemmon2000_air_dAr_ddelta_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d2Ar_ddelta2_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d2Ar_ddelta2, lemmon2000_air_d2Ar_ddelta2_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d3Ar_ddelta3_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d3Ar_ddelta3, lemmon2000_air_d3Ar_ddelta3_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d4Ar_ddelta4_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d4Ar_ddelta4, lemmon2000_air_d4Ar_ddelta4_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8


@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d2Ar_ddeltadtau_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d2Ar_ddeltadtau, lemmon2000_air_d2Ar_ddeltadtau_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8


@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d3Ar_ddeltadtau2_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d3Ar_ddeltadtau2, lemmon2000_air_d3Ar_ddeltadtau2_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d3Ar_ddelta2dtau_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d3Ar_ddelta2dtau, lemmon2000_air_d3Ar_ddelta2dtau_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d4Ar_ddelta2dtau2_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d4Ar_ddelta2dtau2, lemmon2000_air_d4Ar_ddelta2dtau2_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8

@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d4Ar_ddeltadtau3_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d4Ar_ddeltadtau3, lemmon2000_air_d4Ar_ddeltadtau3_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8


@pytest.mark.slow
@pytest.mark.fuzz
def test_lemmon2000_air_d4Ar_ddelta3dtau_vs_naive():
    AARD, std, max_err = func_vs_naive_tester(lemmon2000_air_d4Ar_ddelta3dtau, lemmon2000_air_d4Ar_ddelta3dtau_naive, N=100)
    assert AARD < 1e-13
    assert max_err < 1e-8


def test_lemmon2000_air_rho_dew():
    assert_close(lemmon2000_air_rho_dew(120), 2989.303928859551, rtol=1e-13)

def test_lemmon2000_air_rho_bubble():
    assert_close(lemmon2000_air_rho_bubble(120), 21589.77853554958, rtol=1e-13)


def test_lemmon2000_air_P_dew():
    assert_close(lemmon2000_air_P_dew(100), 567424.1338937179, rtol=1e-14)

def test_lemmon2000_air_P_bubble():
    assert_close(lemmon2000_air_P_bubble(100), 663128.5894402424, rtol=1e-14)


def test_lemmon2000_rho():
    assert_close(lemmon2000_rho(300.0, 1e5), 40.10292351061862, rtol=1e-13)


def test_lemmon2000_P():
    assert_close(lemmon2000_P(300.0, 40.10292351061862), 1e5, rtol=1e-14)


def test_TEOS10_CAAW_derivatives():
    assert_close1d(TEOS10_CAAW_derivatives(200.0)[:3], (1.05493575e-09, -1.525350000000001e-12, -1.13436375e-13), rtol=1e-13)
    assert_close(derivative(lambda T: TEOS10_CAAW_derivatives(T)[-2], 200.0, dx=200*1e-7), TEOS10_CAAW_derivatives(200)[-1], rtol=1e-8)
    assert_close1d(TEOS10_CAAW_derivatives(300.0), (8.019777407407409e-10, -1.9610345679012353e-12, 1.700556378600824e-14, -1.0129827160493834e-16), rtol=1e-13)

def test_TEOS10_CAWW_derivatives():
    assert_close(TEOS10_CAWW_derivatives(200.0)[0], -0.349872634E-5, rtol=1e-9)
    assert_close(TEOS10_CAWW_derivatives(200.0)[1], 0.188025052E-6, rtol=2e-9)
    assert_close(TEOS10_CAWW_derivatives(200.0)[2], -0.124996856E-7, rtol=2e-9)
    assert_close(derivative(lambda T: TEOS10_CAWW_derivatives(T)[-2], 200.0, dx=200*1e-7), TEOS10_CAWW_derivatives(200)[-1], rtol=1e-8)

def test_TEOS10_BAW_derivatives():
    assert_close(TEOS10_BAW_derivatives(200)[0], -0.784874278E-4)
    assert_close(TEOS10_BAW_derivatives(200)[1], 0.848076624E-6)
    assert_close(TEOS10_BAW_derivatives(200)[2], -0.122622146E-7)
    assert_close(derivative(lambda T: TEOS10_BAW_derivatives(T)[-2], 200.0, dx=200*1e-7), TEOS10_BAW_derivatives(200)[-1], rtol=1e-8)

def test_iapws04_Henry_air():
    assert_close(iapws04_Henry_air(300.0), 1.3616423498770563e-10, rtol=1e-13)

def test_iapws04_dHenry_air_dT():
    assert_close(iapws04_dHenry_air_dT(300.0)[0], -1.8759741871455597e-12, rtol=1e-8)
    assert_close(iapws04_dHenry_air_dT(300.0)[1], 1.3616423498770563e-10, rtol=1e-12)
