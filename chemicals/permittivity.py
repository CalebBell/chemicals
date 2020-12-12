# -*- coding: utf-8 -*-
r"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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

This module contains various permittivity calculation routines and dataframes
of coefficients for correlations.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Correlations for Specific Substances
------------------------------------
.. autofunction:: chemicals.permittivity.permittivity_IAPWS
.. autofunction:: chemicals.permittivity.permittivity_CRC

Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: permittivity_data_CRC

    Data from [1]_ with coefficients fit to a polynomial in terms of temperature
    in K. :math:`\epsilon_r = A + BT + CT^2 + DT^3` is the equation, although
    some chemcials only have a constant value.

.. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
   Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.permittivity.permittivity_data_CRC

"""

from __future__ import division

__all__ = ['permittivity_IAPWS', 'permittivity_CRC']

from fluids.numerics import numpy as np
from chemicals.utils import sqrt, PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import register_df_source, data_source

folder = os_path_join(source_path, 'Electrolytes')
register_df_source(folder, 'Permittivity (Dielectric Constant) of Liquids.tsv')

_permittivity_data_loaded = False
def _load_permittivity_data():
    global _permittivity_data_loaded, permittivity_values_CRC, permittivity_data_CRC
    permittivity_data_CRC = data_source('Permittivity (Dielectric Constant) of Liquids.tsv')
    permittivity_values_CRC = np.array(permittivity_data_CRC.values[:, 1:], dtype=float)

if PY37:
    def __getattr__(name):
        if name in ('permittivity_values_CRC', 'permittivity_data_CRC'):
            _load_permittivity_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_permittivity_data()

def permittivity_CRC(T, a, b, c, d):
    r"""
    Return the relative permittivity (epsilon) of a chemical using
    a polynomical equation as in [1]_.

    Parameters
    ----------
    a,b,c,d : float
        Regressed coefficients.

    Notes
    -----
    The permittivity is given by :math:`\epsilon_r = A + BT + CT^2 + DT^3`

    Examples
    --------
    Calculate the permittivity of 4-Nitroaniline:

    >>> permittivity_CRC(450., 487, -1.5, 0.00129, 0.)
    73.225

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
        Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.

    """
    return a + T*(b + T*(c + d*T))

def permittivity_IAPWS(T, rho):
    r'''Calculate the relative permittivity of pure water as a function of.
    temperature and density. Assumes the 1997 IAPWS [1]_ formulation.

    .. math::
        \epsilon(\rho, T) =\frac{1 + A + 5B + (9 + 2A + 18B + A^2 + 10AB +
        9B^2)^{0.5}}{4(1-B)}

    .. math::
        A(\rho, T) = \frac{N_A\mu^2\rho g}{M\epsilon_0 kT}

    .. math::
        B(\rho) = \frac{N_A\alpha\rho}{3M\epsilon_0}

    .. math::
        g(\delta,\tau) = 1 + \sum_{i=1}^{11}n_i\delta^{I_i}\tau^{J_i}
        + n_{12}\delta\left(\frac{647.096}{228}\tau^{-1} - 1\right)^{-1.2}

    .. math::
        \delta = \rho/(322 \text{ kg/m}^3)

    .. math::
        \tau = T/647.096\text{K}

    Parameters
    ----------
    T : float
        Temperature of water [K]
    rho : float
        Mass density of water at T and P [kg/m^3]

    Returns
    -------
    epsilon : float
        Relative permittivity of water at T and rho, [-]

    Notes
    -----
    Validity:

    273.15 < T < 323.15 K for 0 < P < iceVI melting pressure at T or 1000 MPa,
    whichever is smaller.

    323.15 < T < 873.15 K 0 < p < 600 MPa.

    Coefficients and constants (they are optimized away in the function itself):

    ih = [1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10]

    jh = [0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10]

    Nh = [0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396,
          -0.298217036956, -0.108863472196, 0.949327488264E-1,
          -.980469816509E-2, 0.165167634970E-4, 0.937359795772E-4,
          -0.12317921872E-9]

    polarizability = 1.636E-40

    dipole = 6.138E-30

    Examples
    --------
    >>> permittivity_IAPWS(373., 958.46)
    55.565841872697234

    >>> permittivity_IAPWS(650., 40.31090)
    1.2659205723606064

    References
    ----------
    .. [1] IAPWS. 1997. Release on the Static Dielectric Constant of Ordinary
       Water Substance for Temperatures from 238 K to 873 K and Pressures up
       to 1000 MPa.
    '''
#    k = 1.38064852e-23
    # actual molecular dipole moment of water, in C*m
#    dipole2 = 3.7675044000000003e-59 # dipole = 6.138E-30 # but we square it
#    polarizability = 1.636E-40 # actual mean molecular polarizability of water, C^2/J*m^2
#    MW = 0.018015268 # molecular weight of water, kg/mol
#    N_A = 6.0221367e23
    delta = rho*0.003105590062111801 # 1/322.0
#    delta = rho/322.
    T_inv = 1.0/T
    tau = 647.096*T_inv

    g = 1.0  + 0.196096504426E-2*delta*(T*0.0043859649122807015 - 1.0)**-1.2 # 0.00438.. == 1/228.0
#    g = 1.0  + 0.196096504426E-2*delta*(T/228. - 1.0)**-1.2
#    ih = [1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10]
#    jh = [0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10]
#    Nh = [0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396,
#          -0.298217036956, -0.108863472196, 0.949327488264E-1,
#          -.980469816509E-2, 0.165167634970E-4, 0.937359795772E-4,
#          -0.12317921872E-9]
#
#    for h in range(11):
#        g += Nh[h]*delta**ih[h]*tau**jh[h]
    tau_rt = sqrt(tau)
    tau_15 = tau*tau_rt
    tau_25 = tau_15*tau
    tau_2 = tau*tau
    tau_5 = tau_25*tau_25

    g += delta*(-delta*(delta*(-delta*(delta*(delta*(delta*(-1.2317921872000001181e-10
                                                            *delta*delta*delta*tau_5*tau_5
                    + 0.000093735979577200004786*tau_rt)
        + 0.000016516763497000000026*tau_5) - 0.0098046981650899995425*tau_2)
        + 0.094932748826400001341*tau_2) + 0.10886347219599999681*tau_25 + 0.29821703695599999229*tau_15)
        - 0.71469224439599998711*tau_15) + 0.97822448682600005032*sqrt(tau_rt)
        - 0.9577713793750000093*tau + 0.23751179414799999945*tau_25)

    A = rho*g*T_inv*10.302249930663786#N_A*dipole2/(epsilon_0*k*MW)
    B = rho*0.00020588442304500529#Na*polarizability/1(3.*epsilon_0*MW) simplifies to the constant
    epsilon = (1. + A + 5.*B + sqrt(9. + A*(2.0 + A) + B*(18. + 10.*A + 9.*B)
        ))/(4. - 4.*B)
    return epsilon

