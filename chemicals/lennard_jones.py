# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.

This module contains lookup functions and estimation methods for the
parameters molecular diameter `sigma` and the Stockmayer parameter `epsilon`.
These are used for diffusivity calculations. It also contains several
methods for computing the collision integral, another parameter used in the
Lennard-Jones model.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Stockmayer Parameter
--------------------
.. autofunction:: chemicals.lennard_jones.Stockmayer
.. autofunction:: chemicals.lennard_jones.Stockmayer_methods
.. autodata:: chemicals.lennard_jones.Stockmayer_all_methods

Stockmayer Parameter Correlations
---------------------------------
.. autofunction:: chemicals.lennard_jones.epsilon_Flynn
.. autofunction:: chemicals.lennard_jones.epsilon_Bird_Stewart_Lightfoot_critical
.. autofunction:: chemicals.lennard_jones.epsilon_Bird_Stewart_Lightfoot_boiling
.. autofunction:: chemicals.lennard_jones.epsilon_Bird_Stewart_Lightfoot_melting
.. autofunction:: chemicals.lennard_jones.epsilon_Stiel_Thodos
.. autofunction:: chemicals.lennard_jones.epsilon_Tee_Gotoh_Steward_1
.. autofunction:: chemicals.lennard_jones.epsilon_Tee_Gotoh_Steward_2

Molecular Diameter
------------------
.. autofunction:: chemicals.lennard_jones.molecular_diameter
.. autofunction:: chemicals.lennard_jones.molecular_diameter_methods
.. autodata:: chemicals.lennard_jones.molecular_diameter_all_methods

Molecular Diameter Correlations
-------------------------------
.. autofunction:: chemicals.lennard_jones.sigma_Flynn
.. autofunction:: chemicals.lennard_jones.sigma_Bird_Stewart_Lightfoot_critical_2
.. autofunction:: chemicals.lennard_jones.sigma_Bird_Stewart_Lightfoot_critical_1
.. autofunction:: chemicals.lennard_jones.sigma_Bird_Stewart_Lightfoot_boiling
.. autofunction:: chemicals.lennard_jones.sigma_Bird_Stewart_Lightfoot_melting
.. autofunction:: chemicals.lennard_jones.sigma_Stiel_Thodos
.. autofunction:: chemicals.lennard_jones.sigma_Tee_Gotoh_Steward_1
.. autofunction:: chemicals.lennard_jones.sigma_Tee_Gotoh_Steward_2
.. autofunction:: chemicals.lennard_jones.sigma_Silva_Liu_Macedo

Utility Functions
-----------------
.. autofunction:: chemicals.lennard_jones.T_star

"""



__all__ = ['Stockmayer_all_methods', 'Stockmayer_methods', 'Stockmayer',
           'molecular_diameter_all_methods', 'molecular_diameter', 'molecular_diameter_methods',

           'sigma_Flynn',
           'sigma_Bird_Stewart_Lightfoot_critical_2',
           'sigma_Bird_Stewart_Lightfoot_critical_1',
           'sigma_Bird_Stewart_Lightfoot_boiling',
           'sigma_Bird_Stewart_Lightfoot_melting',
           'sigma_Stiel_Thodos', 'sigma_Tee_Gotoh_Steward_1',
           'sigma_Tee_Gotoh_Steward_2',
           'sigma_Silva_Liu_Macedo', 'epsilon_Flynn',
           'epsilon_Bird_Stewart_Lightfoot_critical',
           'epsilon_Bird_Stewart_Lightfoot_boiling',
           'epsilon_Bird_Stewart_Lightfoot_melting', 'epsilon_Stiel_Thodos',
           'epsilon_Tee_Gotoh_Steward_1', 'epsilon_Tee_Gotoh_Steward_2',
           'collision_integral_Neufeld_Janzen_Aziz',
           'As_collision', 'Bs_collision', 'Cs_collision',
           'collision_integral_Kim_Monroe', 'T_star']

import os
from fluids.constants import k
from chemicals.utils import exp, log, sin
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)
#  Register data sources and lazy load them

folder = os_path_join(source_path, 'Viscosity')
register_df_source(folder, 'MagalhaesLJ.tsv')

FLYNN = 'Flynn (1960)'
STIELTHODOS = 'Stiel and Thodos Tc, Zc (1962)'
MAGALHAES = 'Magalhães, Lito, Da Silva, and Silva (2013)'
TEEGOTOSTEWARD1 = 'Tee, Gotoh, and Stewart CSP with Tc (1966)'
TEEGOTOSTEWARD2 = 'Tee, Gotoh, and Stewart CSP with Tc, omega (1966)'
BSLC = 'Bird, Stewart, and Light (2002) critical relation'
BSLB = 'Bird, Stewart, and Light (2002) boiling relation'
BSLM = 'Bird, Stewart, and Light (2002) melting relation'

_LJ_data_loaded = False
def _load_LJ_data():
    global _LJ_data_loaded, LJ_data_Magalhaes, LJ_sources
    LJ_data_Magalhaes = data_source('MagalhaesLJ.tsv')
    _LJ_data_loaded = True
    LJ_sources = {
        MAGALHAES: LJ_data_Magalhaes,
    }

if PY37:
    def __getattr__(name):
        if name in ('LJ_data_Magalhaes', 'LJ_sources'):
            _load_LJ_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_LJ_data()

Stockmayer_all_methods = (MAGALHAES, TEEGOTOSTEWARD2, STIELTHODOS, FLYNN, BSLC,
                          TEEGOTOSTEWARD1, BSLB, BSLM)
'''Tuple of method name keys. See the `Stockmayer` for the actual references'''

def Stockmayer_methods(CASRN=None, Tm=None, Tb=None, Tc=None, Zc=None, omega=None):
    """Return all methods available to obtain the Stockmayer parameter for the
    desired chemical.

    Parameters
    ----------
    CASRN : str, optional
        CASRN [-]
    Tm : float, optional
        Melting temperature of compound [K]
    Tb : float, optional
        Boiling temperature of compound [K]
    Tc : float, optional
        Critical temperature of compound, [K]
    Zc : float, optional
        Critical compressibility of compound, [-]
    omega : float, optional
        Acentric factor of compound, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain `Stockmayer` with the given inputs.

    See Also
    --------
    Stockmayer
    """
    if not _LJ_data_loaded: _load_LJ_data()
    methods = list_available_methods_from_df_dict(LJ_sources, CASRN, 'epsilon')
    if Tc:
        if omega: methods.append(TEEGOTOSTEWARD2)
        if Zc: methods.append(STIELTHODOS)
        methods.append(FLYNN)
        methods.append(BSLC)
        methods.append(TEEGOTOSTEWARD1)
    if Tb: methods.append(BSLB)
    if Tm: methods.append(BSLM)
    return methods

def Stockmayer(CASRN='', Tm=None, Tb=None, Tc=None, Zc=None, omega=None,
               method=None):
    r'''This function handles the retrieval or calculation a chemical's
    Stockmayer parameter. Values are available from one source with lookup
    based on CASRNs, or can be estimated from 7 CSP methods.
    Will automatically select a data source to use if no method is provided;
    returns None if the data is not available.

    Preferred sources are 'Magalhães, Lito, Da Silva, and Silva (2013)' for
    common chemicals which had valies listed in that source, and the CSP method
    `Tee, Gotoh, and Stewart CSP with Tc, omega (1966)` for chemicals which
    don't.

    Examples
    --------
    >>> Stockmayer(CASRN='64-17-5')
    1291.41

    Parameters
    ----------
    CASRN : str, optional
        CASRN [-]
    Tm : float, optional
        Melting temperature of compound [K]
    Tb : float, optional
        Boiling temperature of compound [K]
    Tc : float, optional
        Critical temperature of compound, [K]
    Zc : float, optional
        Critical compressibility of compound, [-]
    omega : float, optional
        Acentric factor of compound, [-]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        Stockmayer_all_methods

    Notes
    -----
    These values are somewhat rough, as they attempt to pigeonhole a chemical
    into L-J behavior.

    The tabulated data is from [2]_, for 322 chemicals.

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    .. [2] Magalhães, Ana L., Patrícia F. Lito, Francisco A. Da Silva, and
       Carlos M. Silva. "Simple and Accurate Correlations for Diffusion
       Coefficients of Solutes in Liquids and Supercritical Fluids over Wide
       Ranges of Temperature and Density." The Journal of Supercritical Fluids
       76 (April 2013): 94-114. doi:10.1016/j.supflu.2013.02.002.
    '''
    if not _LJ_data_loaded: _load_LJ_data()
    if method is not None:
        if method == FLYNN:
            return epsilon_Flynn(Tc)
        elif method == BSLC:
            return epsilon_Bird_Stewart_Lightfoot_critical(Tc)
        elif method == BSLB:
            return epsilon_Bird_Stewart_Lightfoot_boiling(Tb)
        elif method == BSLM:
            return epsilon_Bird_Stewart_Lightfoot_melting(Tm)
        elif method == STIELTHODOS:
            return epsilon_Stiel_Thodos(Tc, Zc)
        elif method == TEEGOTOSTEWARD1:
            return epsilon_Tee_Gotoh_Steward_1(Tc)
        elif method == TEEGOTOSTEWARD2:
            return epsilon_Tee_Gotoh_Steward_2(Tc, omega)
        else:
            return retrieve_from_df_dict(LJ_sources, CASRN, 'epsilon', method)
    else:
        epsilon = retrieve_any_from_df_dict(LJ_sources, CASRN, 'epsilon')
        if epsilon is not None: return epsilon
        if Tc:
            if omega: return epsilon_Tee_Gotoh_Steward_2(Tc, omega)
            if Zc: return epsilon_Stiel_Thodos(Tc, Zc)
            return epsilon_Flynn(Tc)
        if Tb: return epsilon_Bird_Stewart_Lightfoot_boiling(Tb)
        if Tm: return epsilon_Bird_Stewart_Lightfoot_melting(Tm)

TEEGOTOSTEWARD3 = 'Tee, Gotoh, and Stewart CSP with Tc, Pc (1966)'
TEEGOTOSTEWARD4 = 'Tee, Gotoh, and Stewart CSP with Tc, Pc, omega (1966)'
BSLC1 = 'Bird, Stewart, and Light (2002) critical relation with Vc'
BSLC2 = 'Bird, Stewart, and Light (2002) critical relation with Tc, Pc'
STIELTHODOSMD = 'Stiel and Thodos Vc, Zc (1962)'
SILVALIUMACEDO = 'Silva, Liu, and Macedo (1998) critical relation with Tc, Pc'
molecular_diameter_all_methods = (MAGALHAES, TEEGOTOSTEWARD4, SILVALIUMACEDO,
                                  BSLC2, TEEGOTOSTEWARD3, STIELTHODOSMD, FLYNN,
                                  BSLC1, BSLB, BSLM)
'''Tuple of method name keys. See the `molecular_diameter` for the actual references'''

def molecular_diameter_methods(CASRN=None, Tc=None, Pc=None, Vc=None, Zc=None,
                               omega=None, Vm=None, Vb=None):
    """Return all methods available to obtain the molecular diameter for the
    desired chemical.

    Parameters
    ----------
    CASRN : str, optional
        CASRN [-]
    Tc : float, optional
        Critical temperature, [K]
    Pc : float, optional
        Critical pressure, [Pa]
    Vc : float, optional
        Critical volume, [m^3/mol]
    Zc : float, optional
        Critical compressibility, [-]
    omega : float, optional
        Acentric factor of compound, [-]
    Vm : float, optional
        Molar volume of liquid at the melting point of the fluid [K]
    Vb : float, optional
        Molar volume of liquid at the boiling point of the fluid [K]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain `molecular_diameter` with the given inputs.

    See Also
    --------
    molecular_diameter
    """
    methods = list_available_methods_from_df_dict(LJ_sources, CASRN, 'sigma')
    if Tc:
        if Pc:
            if omega: methods.append(TEEGOTOSTEWARD4)
            methods.append(SILVALIUMACEDO)
            methods.append(BSLC2)
            methods.append(TEEGOTOSTEWARD3)
    if Vc:
        if Zc: methods.append(STIELTHODOSMD)
        methods.append(FLYNN)
        methods.append(BSLC1)
    if Vb: methods.append(BSLB)
    if Vm: methods.append(BSLM)
    return methods

def molecular_diameter(CASRN=None, Tc=None, Pc=None, Vc=None, Zc=None, omega=None,
          Vm=None, Vb=None, method=None):
    r'''This function handles the retrieval or calculation a chemical's
    L-J molecular diameter. Values are available from one source with lookup
    based on CASRNs, or can be estimated from 9 CSP methods.
    Will automatically select a data source to use if no method is provided;
    returns None if the data is not available.

    Preferred sources are 'Magalhães, Lito, Da Silva, and Silva (2013)' for
    common chemicals which had valies listed in that source, and the CSP method
    `Tee, Gotoh, and Stewart CSP with Tc, Pc, omega (1966)` for chemicals which
    don't.

    Examples
    --------
    >>> molecular_diameter(CASRN='64-17-5')
    4.23738

    Parameters
    ----------
    CASRN : str, optional
        CASRN [-]
    Tc : float, optional
        Critical temperature, [K]
    Pc : float, optional
        Critical pressure, [Pa]
    Vc : float, optional
        Critical volume, [m^3/mol]
    Zc : float, optional
        Critical compressibility, [-]
    omega : float, optional
        Acentric factor of compound, [-]
    Vm : float, optional
        Molar volume of liquid at the melting point of the fluid [K]
    Vb : float, optional
        Molar volume of liquid at the boiling point of the fluid [K]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        molecular_diameter_all_methods

    Notes
    -----
    These values are somewhat rough, as they attempt to pigeonhole a chemical
    into L-J behavior.

    The tabulated data is from [2]_, for 322 chemicals.

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    .. [2] Magalhães, Ana L., Patrícia F. Lito, Francisco A. Da Silva, and
       Carlos M. Silva. "Simple and Accurate Correlations for Diffusion
       Coefficients of Solutes in Liquids and Supercritical Fluids over Wide
       Ranges of Temperature and Density." The Journal of Supercritical Fluids
       76 (April 2013): 94-114. doi:10.1016/j.supflu.2013.02.002.
    '''
    if not _LJ_data_loaded: _load_LJ_data()
    if method is not None:
        if method == FLYNN:
            return sigma_Flynn(Vc)
        elif method == BSLC1:
            return sigma_Bird_Stewart_Lightfoot_critical_1(Vc)
        elif method == BSLC2:
            return sigma_Bird_Stewart_Lightfoot_critical_2(Tc, Pc)
        elif method == TEEGOTOSTEWARD3:
            return sigma_Tee_Gotoh_Steward_1(Tc, Pc)
        elif method == SILVALIUMACEDO:
            return sigma_Silva_Liu_Macedo(Tc, Pc)
        elif method == BSLB:
            return sigma_Bird_Stewart_Lightfoot_boiling(Vb)
        elif method == BSLM:
            return sigma_Bird_Stewart_Lightfoot_melting(Vm)
        elif method == STIELTHODOSMD:
            return sigma_Stiel_Thodos(Vc, Zc)
        elif method == TEEGOTOSTEWARD4:
            return sigma_Tee_Gotoh_Steward_2(Tc, Pc, omega)
        else:
            return retrieve_from_df_dict(LJ_sources, CASRN, 'sigma', method)
    else:
        epsilon = retrieve_any_from_df_dict(LJ_sources, CASRN, 'sigma')
        if epsilon is not None: return epsilon
        if Tc:
            if Pc:
                if omega: return sigma_Tee_Gotoh_Steward_2(Tc, Pc, omega)
                return sigma_Silva_Liu_Macedo(Tc, Pc)
        if Vc:
            if Zc: return sigma_Stiel_Thodos(Vc, Zc)
            return sigma_Flynn(Vc)
        if Vb: return sigma_Bird_Stewart_Lightfoot_boiling(Vb)
        if Vm: return sigma_Bird_Stewart_Lightfoot_melting(Vm)


### Sigma Lennard-Jones

def sigma_Flynn(Vc):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical volume. CSP method by [1]_ as reported in [2]_.

    .. math::
        \sigma = 0.561(V_c^{1/3})^{5/4}

    Parameters
    ----------
    Vc : float
        Critical volume of fluid [m^3/mol]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Vc is originally in units of mL/mol.

    Examples
    --------
    >>> sigma_Flynn(0.000268)
    5.2506948422196285

    References
    ----------
    .. [1] Flynn, L.W., M.S. thesis, Northwestern Univ., Evanston, Ill. (1960).
    .. [2] Stiel, L. I., and George Thodos. "Lennard-Jones Force Constants
       Predicted from Critical Properties." Journal of Chemical & Engineering
       Data 7, no. 2 (April 1, 1962): 234-36. doi:10.1021/je60013a023
    '''
    Vc = Vc*1E6  # m^3/mol to cm^3/mol
    sigma = 0.561*(Vc**(1/3.))**1.2
    return sigma

def sigma_Bird_Stewart_Lightfoot_critical_2(Tc, Pc):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical temperature and pressure. CSP method by [1]_.

    .. math::
        \sigma = 2.44(T_c/P_c)^{1/3}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Original units of critical pressure are atmospheres.

    Examples
    --------
    >>> sigma_Bird_Stewart_Lightfoot_critical_2(560.1, 4550000)
    5.658657684653222

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    Pc = Pc/101325.
    sigma = 2.44*(Tc/Pc)**(1/3.0)
    return sigma

def sigma_Bird_Stewart_Lightfoot_critical_1(Vc):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical volume. CSP method by [1]_.

    .. math::
        \sigma = 0.841 V_c^{1/3}

    Parameters
    ----------
    Vc : float
        Critical volume of fluid [m^3/mol]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Original units of Vc are mL/mol.

    Examples
    --------
    >>> sigma_Bird_Stewart_Lightfoot_critical_1(0.000268)
    5.422184116631474

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    Vc = Vc*1E6  # m^3/mol to cm^3/mol
    sigma = 0.841*Vc**(1/3.0)
    return sigma

def sigma_Bird_Stewart_Lightfoot_boiling(Vb):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses molar volume of liquid at boiling. CSP method by [1]_.

    .. math::
        \sigma = 1.166V_{b,liq}^{1/3}

    Parameters
    ----------
    Vb : float
        Boiling molar volume of liquid [m^3/mol]

    Returns
    -------
    sigma : float
        Lennard-Jones collision integral, [Angstrom]

    Notes
    -----
    Original units of Vb are mL/mol.

    Examples
    --------
    >>> sigma_Bird_Stewart_Lightfoot_boiling(0.0001015)
    5.439018856944655

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    Vb = Vb*1E6
    sigma = 1.166*Vb**(1/3.0)
    return sigma

def sigma_Bird_Stewart_Lightfoot_melting(Vm):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses molar volume of a liquid at its melting point. CSP method by [1]_.

    .. math::
        \sigma = 1.222 V_{m,sol}^{1/3}

    Parameters
    ----------
    Vm : float
        Melting molar volume of a liquid at its melting point [m^3/mol]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Original units of Vm are mL/mol.

    Examples
    --------
    >>> sigma_Bird_Stewart_Lightfoot_melting(8.8e-05)
    5.435407341351406

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    Vm = Vm*1E6
    sigma = 1.222*Vm**(1/3.)
    return sigma

def sigma_Stiel_Thodos(Vc, Zc):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical volume and compressibility. CSP method by [1]_.

    .. math::
        \sigma = 0.1866 V_c^{1/3} Z_c^{-6/5}

    Parameters
    ----------
    Vc : float
        Critical volume of fluid [m^3/mol]
    Zc : float
        Critical compressibility of fluid, [-]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Vc is originally in units of mL/mol.

    Examples
    --------
    Monofluorobenzene

    >>> sigma_Stiel_Thodos(0.000271, 0.265)
    5.94300853971033

    References
    ----------
    .. [1] Stiel, L. I., and George Thodos. "Lennard-Jones Force Constants
       Predicted from Critical Properties." Journal of Chemical & Engineering
       Data 7, no. 2 (April 1, 1962): 234-36. doi:10.1021/je60013a023
    '''
    Vc = Vc*1E6
    sigma = 0.1866*Vc**(1/3.0)*Zc**(-1.2)
    return sigma


def sigma_Tee_Gotoh_Steward_1(Tc, Pc):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical temperature and pressure. CSP method by [1]_.

    .. math::
        \sigma = 2.3647 \left(\frac{T_c}{P_c}\right)^{1/3}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Original units of Pc are atm. Further regressions with other parameters
    were performed in [1]_ but are not included here, except for
    `sigma_Tee_Gotoh_Steward_2`.

    Examples
    --------
    >>> sigma_Tee_Gotoh_Steward_1(560.1, 4550000)
    5.48402779790962

    References
    ----------
    .. [1] Tee, L. S., Sukehiro Gotoh, and W. E. Stewart. "Molecular
       Parameters for Normal Fluids. Lennard-Jones 12-6 Potential." Industrial
       & Engineering Chemistry Fundamentals 5, no. 3 (August 1, 1966): 356-63.
       doi:10.1021/i160019a011
    '''
    Pc = Pc/101325.
    sigma = 2.3647*(Tc/Pc)**(1/3.)
    return sigma

def sigma_Tee_Gotoh_Steward_2(Tc, Pc, omega):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical temperature, pressure, and acentric factor. CSP method by
    [1]_.

    .. math::
        \sigma = (2.3551 - 0.0874\omega)\left(\frac{T_c}{P_c}\right)^{1/3}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Original units of Pc are atm. Further regressions with other parameters
    were performed in [1]_ but are not included here, except for
    `sigma_Tee_Gotoh_Steward_1`.

    Examples
    --------
    >>> sigma_Tee_Gotoh_Steward_2(560.1, 4550000, 0.245)
    5.412104867264477

    References
    ----------
    .. [1] Tee, L. S., Sukehiro Gotoh, and W. E. Stewart. "Molecular Parameters
       for Normal Fluids. Lennard-Jones 12-6 Potential." Industrial
       & Engineering Chemistry Fundamentals 5, no. 3 (August 1, 1966): 356-63.
       doi:10.1021/i160019a011
    '''
    Pc = Pc/101325.
    sigma = (2.3551-0.0874*omega)*(Tc/Pc)**(1/3.)
    return sigma

def sigma_Silva_Liu_Macedo(Tc, Pc):
    r'''Calculates Lennard-Jones molecular diameter.
    Uses critical temperature and pressure. CSP method by [1]_.

    .. math::
        \sigma_{LJ}^3 = 0.17791 + 11.779 \left( \frac{T_c}{P_c}\right)
        - 0.049029\left( \frac{T_c}{P_c}\right)^2

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    sigma : float
        Lennard-Jones molecular diameter, [Angstrom]

    Notes
    -----
    Pc is originally in bar. An excellent paper. None is
    returned if the polynomial returns a negative number, as in the case of
    1029.13 K and 3.83 bar.

    Examples
    --------
    >>> sigma_Silva_Liu_Macedo(560.1, 4550000)
    5.164483998730177

    References
    ----------
    .. [1] Silva, Carlos M., Hongqin Liu, and Eugenia A. Macedo. "Models for
       Self-Diffusion Coefficients of Dense Fluids, Including Hydrogen-Bonding
       Substances." Chemical Engineering Science 53, no. 13 (July 1, 1998):
       2423-29. doi:10.1016/S0009-2509(98)00037-2
    '''
    Pc = Pc/1E5  # Pa to bar
    term = 0.17791 + 11.779*(Tc/Pc) - 0.049029 * (Tc/Pc)**2
    if term < 0:
        sigma = None
    else:
        sigma = (term)**(1/3.)
    return sigma

### epsilon Lennard-Jones

def epsilon_Flynn(Tc):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses critical temperature. CSP method by [1]_ as reported in [2]_.

    .. math::
        \epsilon/k = 1.77 T_c^{5/6}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----

    Examples
    --------
    >>> epsilon_Flynn(560.1)
    345.2984087011443

    References
    ----------
    .. [1] Flynn, L.W., M.S. thesis, Northwestern Univ., Evanston, Ill. (1960).
    .. [2] Stiel, L. I., and George Thodos. "Lennard-Jones Force Constants
       Predicted from Critical Properties." Journal of Chemical & Engineering
       Data 7, no. 2 (April 1, 1962): 234-36. doi:10.1021/je60013a023
    '''
    epsilon_k = 1.77*Tc**(5/6.)
    return epsilon_k

def epsilon_Bird_Stewart_Lightfoot_critical(Tc):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses critical temperature. CSP method by [1]_.

    .. math::
        \epsilon/k = 0.77T_c

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----

    Examples
    --------
    >>> epsilon_Bird_Stewart_Lightfoot_critical(560.1)
    431.27700000000004

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    epsilon_k = 0.77*Tc
    return epsilon_k

def epsilon_Bird_Stewart_Lightfoot_boiling(Tb):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses boiling temperature. CSP method by [1]_.

    .. math::
        \epsilon/k = 1.15 T_b

    Parameters
    ----------
    Tb : float
        Boiling temperature [K]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----

    Examples
    --------
    >>> epsilon_Bird_Stewart_Lightfoot_boiling(357.85)
    411.5275

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    epsilon_k = 1.15*Tb
    return epsilon_k

def epsilon_Bird_Stewart_Lightfoot_melting(Tm):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses melting temperature. CSP method by [1]_.

    .. math::
        \epsilon/k = 1.92T_m

    Parameters
    ----------
    Tm : float
        Melting temperature [K]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----

    Examples
    --------
    >>> epsilon_Bird_Stewart_Lightfoot_melting(231.15)
    443.808

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    epsilon_k = 1.92*Tm
    return epsilon_k

def epsilon_Stiel_Thodos(Tc, Zc):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses Critical temperature and critical compressibility. CSP method by [1]_.

    .. math::
        \epsilon/k = 65.3 T_c Z_c^{3.6}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Zc : float
        Critical compressibility of fluid, [-]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----

    Examples
    --------
    Fluorobenzene

    >>> epsilon_Stiel_Thodos(358.5, 0.265)
    196.3755830305783

    References
    ----------
    .. [1] Stiel, L. I., and George Thodos. "Lennard-Jones Force Constants
       Predicted from Critical Properties." Journal of Chemical & Engineering
       Data 7, no. 2 (April 1, 1962): 234-36. doi:10.1021/je60013a023
    '''
    epsilon_k = 65.3*Tc*Zc**3.6
    return epsilon_k

def epsilon_Tee_Gotoh_Steward_1(Tc):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses Critical temperature. CSP method by [1]_.

    .. math::
        \epsilon/k = 0.7740T_c

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----
    Further regressions with other parameters were performed in [1]_ but are
    not included here, except for `epsilon_Tee_Gotoh_Steward_2`.

    Examples
    --------
    >>> epsilon_Tee_Gotoh_Steward_1(560.1)
    433.5174

    References
    ----------
    .. [1] Tee, L. S., Sukehiro Gotoh, and W. E. Stewart. "Molecular Parameters
       for Normal Fluids. Lennard-Jones 12-6 Potential." Industrial &
       Engineering Chemistry Fundamentals 5, no. 3 (August 1, 1966): 356-63.
       doi:10.1021/i160019a011
    '''
    epsilon_k = 0.7740*Tc
    return epsilon_k

def epsilon_Tee_Gotoh_Steward_2(Tc, omega):
    r'''Calculates Lennard-Jones depth of potential-energy minimum.
    Uses critical temperature and acentric factor. CSP method by [1]_.

    .. math::
        \epsilon/k = (0.7915 + 0.1693 \omega)T_c

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    epsilon_k : float
        Lennard-Jones depth of potential-energy minimum over k, [K]

    Notes
    -----
    Further regressions with other parameters were performed in [1]_ but are
    not included here, except for `epsilon_Tee_Gotoh_Steward_1`.

    Examples
    --------
    >>> epsilon_Tee_Gotoh_Steward_2(560.1, 0.245)
    466.55125785

    References
    ----------
    .. [1] Tee, L. S., Sukehiro Gotoh, and W. E. Stewart. "Molecular Parameters
       for Normal Fluids. Lennard-Jones 12-6 Potential." Industrial &
       Engineering Chemistry Fundamentals 5, no. 3 (August 1, 1966): 356-63.
       doi:10.1021/i160019a011
    '''
    epsilon_k = (0.7915 + 0.1693*omega)*Tc
    return epsilon_k

### Collision Integral

def collision_integral_Neufeld_Janzen_Aziz(T_star, l=1, s=1):
    r'''Calculates Lennard-Jones collision integral for any of 16 values of
    (l,j) for the wide range of 0.3 < T_star < 100. Values are accurate to
    0.1 % of actual values, but the calculation of actual values is
    computationally intensive and so these simplifications are used, developed
    in [1]_.

    .. math::
        \Omega_D = \frac{A}{T^{*B}} + \frac{C}{\exp(DT^*)} +
        \frac{E}{\exp(FT^{*})} + \frac{G}{\exp(HT^*)} + RT^{*B}\sin(ST^{*W}-P)

    Parameters
    ----------
    Tstar : float
        Reduced temperature of the fluid [-]
    l : int
        term
    s : int
        term

    Returns
    -------
    Omega : float
        Collision integral of A and B

    Notes
    -----
    Acceptable pairs of (l,s) are (1, 1), (1, 2), (1, 3), (1, 4), (1, 5),
    (1, 6), (1, 7), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 3), (3, 4),
    (3, 5), and (4, 4).

    .. math::
        T^* = \frac{k_b T}{\epsilon}

    Results are very similar to those of the more modern formulation,
    `collision_integral_Kim_Monroe`.

    Calculations begin to yield overflow errors in some values of (l, 2) after
    T_star = 75, beginning with (1, 7). Also susceptible are (1, 5) and (1, 6).

    Examples
    --------
    >>> collision_integral_Neufeld_Janzen_Aziz(100, 1, 1)
    0.516717697672334

    References
    ----------
    .. [1] Neufeld, Philip D., A. R. Janzen, and R. A. Aziz. "Empirical
       Equations to Calculate 16 of the Transport Collision Integrals
       Omega(l, S)* for the Lennard-Jones (12-6) Potential." The Journal of
       Chemical Physics 57, no. 3 (August 1, 1972): 1100-1102.
       doi:10.1063/1.1678363
    '''
    if l == 1 and s == 1:
        A, B, C, D, E, F, G, H = 1.06036, 0.1561, 0.193, 0.47635, 1.03587, 1.52996, 1.76474, 3.89411
    elif l == 1 and s == 2:
        A, B, C, D, E, F, G, H = 1.0022, 0.1553, 0.16105, 0.72751, 0.86125, 2.06848, 1.95162, 4.84492
    elif l == 1 and s == 3:
        A, B, C, D, E, F, R, S, W, P = 0.96573, 0.15611, 0.44067, 1.5242, 2.38981, 5.08063, -0.0005373, 19.2866, -1.30775, 6.58711
    elif l == 1 and s == 4:
        A, B, C, D, E, F, R, S, W, P = 0.93447, 0.15578, 0.39478, 1.85761, 2.45988, 6.15727, 0.0004246, 12.988, -1.36399, 3.3329
    elif l == 1 and s == 5:
        A, B, C, D, E, F, R, S, W, P = 0.90972, 0.15565, 0.35967, 2.18528, 2.45169, 7.17936, -0.0003814, 9.38191, 0.14025, 9.93802
    elif l == 1 and s == 6:
        A, B, C, D, E, F, R, S, W, P = 0.88928, 0.15562, 0.33305, 2.51303, 2.36298, 8.1169, -0.0004649, 9.86928, 0.12851, 9.82414
    elif l == 1 and s == 7:
        A, B, C, D, E, F, R, S, W, P = 0.87208, 0.15568, 0.36583, 3.01399, 2.70659, 9.9231, -0.0004902, 10.2274, 0.12306, 9.97712
    elif l == 2 and s == 2:
        A, B, C, D, E, F, R, S, W, P = 1.16145, 0.14874, 0.52487, 0.7732, 2.16178, 2.43787, -0.0006435, 18.0323, -0.7683, 7.27371
    elif l == 2 and s == 3:
        A, B, C, D, E, F, R, S, W, P = 1.11521, 0.14796, 0.44844, 0.99548, 2.30009, 3.06031, 0.0004565, 38.5868, -0.69403, 2.56375
    elif l == 2 and s == 4:
        A, B, C, D, E, F, R, S, W, P = 1.08228, 0.14807, 0.47128, 1.31596, 2.42738, 3.90018, -0.0005623, 3.08449, 0.28271, 3.22871
    elif l == 2 and s == 5:
        A, B, C, D, E, F, R, S, W, P = 1.05581, 0.14822, 0.51203, 1.67007, 2.57317, 4.85939, -0.000712, 4.7121, 0.2173, 4.7353
    elif l == 2 and s == 6:
        A, B, C, D, E, F, R, S, W, P = 1.03358, 0.14834, 0.53928, 2.01942, 2.7235, 5.84817, -0.0008576, 7.66012, 0.15493, 7.6011
    elif l == 3 and s == 3:
        A, B, C, D, E, F, G, H, R, S, W, P = 1.05567, 0.1498, 0.30887, 0.86437, 1.35766, 2.44123, 1.2903, 5.55734, 0.0002339, 57.7757, -1.0898, 6.9475
    elif l == 3 and s == 4:
        A, B, C, D, E, F, R, S, W, P = 1.02621, 0.1505, 0.55381, 1.4007, 2.06176, 4.26234, 0.0005227, 11.3331, -0.8209, 3.87185
    elif l == 3 and s == 5:
        A, B, C, D, E, F, R, S, W, P = 0.99958, 0.15029, 0.50441, 1.64304, 2.06947, 4.87712, -0.0005184, 3.45031, 0.26821, 3.73348
    elif l == 4 and s == 4:
        A, B, C, D, E, F, R, S, W, P = 1.12007, 0.14578, 0.53347, 1.11986, 2.28803, 3.27567, 0.0007427, 21.048, -0.28759, 6.69149
    else:
        raise ValueError('Input values of l and s are not supported')
    omega = A/T_star**B + C*exp(-D*T_star) + E*exp(-F*T_star)
    if (l == 1 and (s == 1 or s == 2)) or (s == 3 and l == 3):
        omega += G*exp(-H*T_star)
    if not (l == 1  and (s == 1 or s == 2)):
        omega += R*T_star**B*sin(S*T_star**W - P)
    return omega

As_collision = {(1, 1): -1.10367290,
                (1, 2): 1.35555540,
                (1, 3): 1.06771150,
                (1, 4): 0.80959899,
                (1, 5): 0.74128322,
                (1, 6): 0.80998324,
                (1, 7): 0.81808091,
                (2, 2): -0.92032979,
                (2, 3): 2.59557990,
                (2, 4): 1.60427450,
                (2, 5): 0.82064641,
                (2, 6): 0.79413652,
                (3, 3): 1.26304910,
                (3, 4): 2.21146360,
                (3, 5): 1.50498090,
                (4, 4): 2.62223930
                }

Bs_collision = {
(1, 1): [2.6431984,0.0060432255,-0.15158773,0.054237938,-0.0090468682,0.0006174200700],
(1, 2): [-0.44668594,0.42734391,-0.16036459,0.031461648,-0.0032587575,0.0001386025700],
(1, 3): [-0.1394539,0.17696362,-0.026252211,-0.0043814141,0.00167521,-0.0001438280100],
(1, 4): [0.1293817,0.059760309,0.0071109469,-0.0063851124,0.0010498938,-0.0000581492570],
(1, 5): [0.1778885,0.027398438,0.0076254248,-0.0031650182,0.0003278652,-0.0000092890016],
(1, 6): [0.073071217,0.034607908,-0.0011457199,0.000281986,-0.0002006054,0.0000214464830],
(1, 7): [0.044232851,0.029750283,-0.0022011682,0.0006326412,-0.0001755553,0.0000142557040],
(2, 2): [2.3508044,0.50110649,-0.47193769,0.15806367,-0.026367184,0.0018120118000],
(2, 3): [-1.8569443,0.96985775,-0.39888526,0.090063692,-0.010918991,0.0005664679700],
(2, 4): [-0.67406115,0.42671907,-0.10177069,0.0006185714,0.0031225358,-0.0003520605100],
(2, 5): [0.23195128,0.12233793,0.013891578,-0.020903423,0.0046715462,-0.0003520430300],
(2, 6): [0.23766123,0.077125802,0.013060901,-0.010982362,0.0018034505,-0.0000959825710],
(3, 3): [-0.36104243,0.68116214,-0.36401583,0.10500196,-0.016400134,0.0010880886000],
(3, 4): [-1.4743107,0.64918549,-0.24075196,0.051820149,-0.0060565396,0.0002981232600],
(3, 5): [-0.64335529,0.3261704,-0.082126072,0.0059682011,0.0010269488,-0.0001595725200],
(4, 4): [-1.9158462,1.016638,-0.43355278,0.10496591,-0.013951104,0.0008004853400]
}

Cs_collision = {
    (1, 1): [1.6690746, -0.6914589, 0.15502132, -0.020642189, 0.001540207700, -0.000049729535],
    (1, 2): [-0.47499422, 0.14482036, -0.032158368, 0.0044357933, -0.00034138118, 0.000011259742],
    (1, 3): [-0.25258689, 0.059709197, -0.013332695, 0.0019619285, -0.000160630760, 0.0000055804557],
    (1, 4): [-0.045055948, -0.022642753, 0.0056672308, -0.0006570876, 0.000040733113, -0.0000010820157],
    (1, 5): [0.0013668724, -0.041730962, 0.010378923, -0.0013492954, 0.000096963599, -0.0000030307552],
    (1, 6): [-0.071180849, -0.012738119, 0.0038582834, -0.0004706043, 0.000030466929, -0.00000085305576],
    (1, 7): [-0.089417548, -0.0051856424, 0.0021882143, -0.0002487447, 0.000013745859, -0.00000030285365],
    (2, 2): [1.6330213, -0.69795156, 0.16096572, -0.02210944, 0.0017031434, -0.000056699986],
    (2, 3): [-1.4586197, 0.52947262, -0.11946363, 0.016264589, -0.0012354315, 0.000040366357],
    (2, 4): [-0.62774499, 0.20700644, -0.04760169, 0.0067153792, -0.00052706167, 0.000017705708],
    (2, 5): [0.039184885, -0.057316906, 0.012794497, -0.0015336449, 0.00010241454, -0.0000029975563],
    (2, 6): [0.050470266, -0.062621672, 0.014326724, -0.0017806541, 0.00012353365, -0.0000037501381],
    (3, 3): [-0.33227158, 0.079723851, -0.015470355, 0.0018686705, -0.00012179945, 0.0000032594587],
    (3, 4): [-1.1942554, 0.43000688, -0.097525871, 0.013399366, -0.0010283777, 0.000033956674],
    (3, 5): [-0.60014514, 0.19764859, -0.045212434, 0.0063650284, -0.00049991689, 0.000016833944],
    (4, 4): [-1.4676253, 0.53048161, -0.11909781, 0.016123847, -0.0012174905, 0.0000395451]
}

def collision_integral_Kim_Monroe(T_star, l=1, s=1):
    r'''Calculates Lennard-Jones collision integral for any of 16 values of
    (l,j) for the wide range of 0.3 < T_star < 400. Values are accurate to
    0.007 % of actual values, but the calculation of actual values is
    computationally intensive and so these simplifications are used, developed
    in [1]_.

    .. math::
        \Omega^{(l,s)*} = A^{(l,s)} + \sum_{k=1}^6 \left[ \frac{B_k^{(l,s)}}
        {(T^*)^k} + C_k^{(l,s)} (\ln T^*)^k \right]

    Parameters
    ----------
    Tstar : float
        Reduced temperature of the fluid [-]
    l : int
        term
    s : int
        term

    Returns
    -------
    Omega : float
        Collision integral of A and B

    Notes
    -----
    Acceptable pairs of (l,s) are (1, 1), (1, 2), (1, 3), (1, 4), (1, 5),
    (1, 6), (1, 7), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 3), (3, 4),
    (3, 5), and (4, 4).

    .. math::
        T^* = \frac{k_b T}{\epsilon}

    Examples
    --------
    >>> collision_integral_Kim_Monroe(400, 1, 1)
    0.4141818082392228

    References
    ----------
    .. [1] Kim, Sun Ung, and Charles W. Monroe. "High-Accuracy Calculations of
       Sixteen Collision Integrals for Lennard-Jones (12-6) Gases and Their
       Interpolation to Parameterize Neon, Argon, and Krypton." Journal of
       Computational Physics 273 (September 15, 2014): 358-73.
       doi:10.1016/j.jcp.2014.05.018.
    '''
    if (l, s) not in As_collision:
        raise ValueError('Input values of l and s are not supported')
    omega = As_collision[(l, s)]
    Bs = Bs_collision[(l, s)]
    Cs = Cs_collision[(l, s)]
    for ki in range(6):
        omega += Bs[ki] / T_star ** (ki + 1) + Cs[ki] * log(T_star) ** (ki + 1)
    return omega

### Misc

def T_star(T, epsilon_k=None, epsilon=None):
    r'''This function calculates the parameter `T_star` as needed in performing
    collision integral calculations.

    .. math::
        T^* = \frac{kT}{\epsilon}

    Examples
    --------
    >>> T_star(T=318.2, epsilon_k=308.43)
    1.0316765554582887

    Parameters
    ----------
    epsilon_k : float, optional
        Lennard-Jones depth of potential-energy minimum over k, [K]
    epsilon : float, optional
        Lennard-Jones depth of potential-energy minimum [J]

    Returns
    -------
    T_star : float
        Dimentionless temperature for calculating collision integral, [-]

    Notes
    -----
    Tabulated values are normally listed as epsilon/k. k is the Boltzman
    constant, with units of J/K.

    References
    ----------
    .. [1] Bird, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot.
       Transport Phenomena, Revised 2nd Edition. New York:
       John Wiley & Sons, Inc., 2006
    '''
    if epsilon_k:
        _Tstar = T/(epsilon_k)
    elif epsilon:
        _Tstar = k*T/epsilon
    else:
        raise ValueError('Either epsilon/k or epsilon must be provided')
    return _Tstar
