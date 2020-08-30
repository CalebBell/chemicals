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

This module contains a lookup function, a definition function, and correlations
for a chemical's `acentric` factor, normally given the variable :math:`\omega`.

A similar variable called the `stiel polar factor` can be calculated from its
definition as well.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Lookup Functions
----------------
.. autofunction:: chemicals.acentric.omega
.. autofunction:: chemicals.acentric.omega_methods
.. autodata:: chemicals.acentric.omega_all_methods

Definitions
-----------
.. autofunction:: chemicals.acentric.omega_definition
.. autofunction:: chemicals.acentric.Stiel_polar_factor

Correlations
------------
.. autofunction:: chemicals.acentric.LK_omega
"""

__all__ = ['omega', 'LK_omega', 'Stiel_polar_factor',
           'omega_methods', 'omega_all_methods', 'omega_definition']

from chemicals.utils import log, log10
from chemicals import critical
from chemicals.data_reader import (retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)

omega_all_methods = ('PSRK', 'PD', 'YAWS')
'''Tuple of method name keys. See the `omega` for the actual references'''

def omega_methods(CASRN):
    """Return all methods available for obtaining omega for the desired
    chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain omega with the given inputs.

    See Also
    --------
    omega
    """
    return list_available_methods_from_df_dict(critical.omega_sources, CASRN, 'omega')

def omega(CASRN, method=None):
    r'''Retrieve a chemical's acentric factor, `omega`.

    Automatically select a method to use if no method is provided;
    returns None if the data is not available.

    .. math::
        \omega \equiv -\log_{10}\left[\lim_{T/T_c=0.7}(P^{sat}/P_c)\right]-1.0

    Examples
    --------
    >>> omega(CASRN='64-17-5')
    0.635

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    omega : float
        Acentric factor of compound

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'PSRK', 'PD', or 'YAWS'.
        All valid values are also held in the variable `omega_all_methods`.

    Notes
    -----
    A total of three sources are available for this function. They are:

        * 'PSRK', a compillation of experimental and estimated data published 
          in the Appendix of [15]_, the fourth revision of the PSRK model.
        * 'PD', an older compillation of
          data published in (Passut & Danner, 1973) [16]_.
        * 'YAWS', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [17]_.

    See Also
    --------
    omega_methods

    References
    ----------
    .. [1] Pitzer, K. S., D. Z. Lippmann, R. F. Curl, C. M. Huggins, and
       D. E. Petersen: The Volumetric and Thermodynamic Properties of Fluids.
       II. Compressibility Factor, Vapor Pressure and Entropy of Vaporization.
       J. Am. Chem. Soc., 77: 3433 (1955).
    .. [2] Horstmann, Sven, Anna Jabłoniec, Jörg Krafczyk, Kai Fischer, and
       Jürgen Gmehling. "PSRK Group Contribution Equation of State:
       Comprehensive Revision and Extension IV, Including Critical Constants
       and Α-Function Parameters for 1000 Components." Fluid Phase Equilibria
       227, no. 2 (January 25, 2005): 157-64. doi:10.1016/j.fluid.2004.11.002.
    .. [3] Passut, Charles A., and Ronald P. Danner. "Acentric Factor. A
       Valuable Correlating Parameter for the Properties of Hydrocarbons."
       Industrial & Engineering Chemistry Process Design and Development 12,
       no. 3 (July 1, 1973): 365-68. doi:10.1021/i260047a026.
    .. [4] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    if method:
        return retrieve_from_df_dict(critical.omega_sources, CASRN, 'omega', method) 
    else:
        return retrieve_any_from_df_dict(critical.omega_sources, CASRN, 'omega') 

def omega_definition(Psat, Pc):
    r'''Returns the acentric factor of a fluid according to its fundamental
    definition using the vapor pressure at a reduced temperature of 0.7Tc.

    .. math::
        \omega \equiv -\log_{10}\left[\lim_{T/T_c=0.7}(P^{sat}/P_c)\right]-1.0

    Parameters
    ----------
    Psat : float
        Vapor pressure of the fluid at a reduced temperature of 0.7 [Pa]
    Pc : float
        Critical pressure of the fluid [Pa]

    Returns
    -------
    omega : float
        Acentric factor of the fluid [-]

    Notes
    -----

    Examples
    --------
    Water

    >>> omega_definition(999542, 22048320.0)
    0.3435744558761711

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return -log10(Psat/Pc) - 1.0

def LK_omega(Tb, Tc, Pc):
    r'''Estimates the acentric factor of a fluid using a correlation in [1]_.

    .. math::
        \omega = \frac{\ln P_{br}^{sat} - 5.92714 + 6.09648/T_{br} + 1.28862
        \ln T_{br} -0.169347T_{br}^6}
        {15.2518 - 15.6875/T_{br} - 13.4721 \ln T_{br} + 0.43577 T_{br}^6}

    Parameters
    ----------
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]

    Returns
    -------
    omega : float
        Acentric factor of the fluid [-]

    Notes
    -----
    The units of the above equation are atmosphere and Kelvin; values are 
    converted internally.

    Examples
    --------
    Isopropylbenzene, from Reid (1987).

    >>> LK_omega(425.6, 631.1, 32.1E5)
    0.32544249926397856

    References
    ----------
    .. [1] Lee, Byung Ik, and Michael G. Kesler. "A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States." AIChE Journal
       21, no. 3 (1975): 510-527. doi:10.1002/aic.690210313.
    '''
    T_br = Tb/Tc
    log_T_br = log(T_br)
    T_br_6 = T_br*T_br
    T_br_6 *= T_br_6*T_br_6
    omega = (log(101325.0/Pc) - 5.92714 + 6.09648/T_br + 1.28862*log_T_br -
             0.169347*T_br_6)/(15.2518 - 15.6875/T_br - 13.4721*log_T_br +
             0.43577*T_br_6)
    return omega

def Stiel_polar_factor(Psat, Pc, omega):
    r'''This function handles the calculation of a chemical's Stiel Polar
    factor, directly through the definition of Stiel-polar factor.
    Requires the vapor pressure `Psat` at a reduced temperature of 0.6,
    the critical pressure `Pc`, and the acentric factor `omega`.

    .. math::
        x = \log P_r|_{T_r=0.6} + 1.70 \omega + 1.552

    Parameters
    ----------
    Psat : float
        Vapor pressure of fluid at a reduced temperature of 0.6 [Pa]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor of the fluid [-]

    Returns
    -------
    factor : float
        Stiel polar factor of compound

    Notes
    -----
    A few points have also been published in [2]_, which may be used for
    comparison. Currently this is only used for a surface tension correlation.

    Examples
    --------
    Calculating the factor for water:
    
    >>> Stiel_polar_factor(Psat=169745, Pc=22048321.0, omega=0.344)
    0.02322146744772713

    References
    ----------
    .. [1] Halm, Roland L., and Leonard I. Stiel. "A Fourth Parameter for the
       Vapor Pressure and Entropy of Vaporization of Polar Fluids." AIChE
       Journal 13, no. 2 (1967): 351-355. doi:10.1002/aic.690130228.
    .. [2] D, Kukoljac Miloš, and Grozdanić Dušan K. "New Values of the
       Polarity Factor." Journal of the Serbian Chemical Society 65, no. 12
       (January 1, 2000). http://www.shd.org.rs/JSCS/Vol65/No12-Pdf/JSCS12-07.pdf
    '''
    Pr = Psat/Pc
    factor = log10(Pr) + 1.70*omega + 1.552
    return factor
