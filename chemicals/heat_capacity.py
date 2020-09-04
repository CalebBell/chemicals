# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell
<Caleb.Andrew.Bell@gmail.com>
Copyright (C) 2020 Yoel Rene Cortes-Pena
<yoelcortes@gmail.com>

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

This module contains many heat capacity model equations, heat capacity estimation
equations, enthalpy and entropy integrals of those heat capacity equations, 
enthalpy/entropy flash initialization routines, and many dataframes of 
coefficients.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:
    
Gas Heat Capacity Model Equations
---------------------------------
.. autofunction:: chemicals.heat_capacity.TRCCp
.. autofunction:: chemicals.heat_capacity.TRCCp_integral
.. autofunction:: chemicals.heat_capacity.TRCCp_integral_over_T

Gas Heat Capacity Estimation Models
-----------------------------------
.. autofunction:: chemicals.heat_capacity.Poling
.. autofunction:: chemicals.heat_capacity.Poling_integral
.. autofunction:: chemicals.heat_capacity.Poling_integral_over_T
.. autofunction:: chemicals.heat_capacity.Lastovka_Shaw
.. autofunction:: chemicals.heat_capacity.Lastovka_Shaw_integral
.. autofunction:: chemicals.heat_capacity.Lastovka_Shaw_integral_over_T
.. autofunction:: chemicals.heat_capacity.Lastovka_Shaw_T_for_Hm
.. autofunction:: chemicals.heat_capacity.Lastovka_Shaw_T_for_Sm
.. autofunction:: chemicals.heat_capacity.Lastovka_Shaw_term_A

Liquid Heat Capacity Model Equations
------------------------------------
.. autofunction:: chemicals.heat_capacity.Zabransky_quasi_polynomial
.. autofunction:: chemicals.heat_capacity.Zabransky_quasi_polynomial_integral
.. autofunction:: chemicals.heat_capacity.Zabransky_quasi_polynomial_integral_over_T
.. autofunction:: chemicals.heat_capacity.Zabransky_cubic
.. autofunction:: chemicals.heat_capacity.Zabransky_cubic_integral
.. autofunction:: chemicals.heat_capacity.Zabransky_cubic_integral_over_T
.. autofunction:: chemicals.heat_capacity.Zabransky_quasi_polynomial
.. autoclass:: chemicals.heat_capacity.ZabranskySpline
.. autoclass:: chemicals.heat_capacity.ZabranskyQuasipolynomial

Liquid Heat Capacity Estimation Models
--------------------------------------
.. autofunction:: chemicals.heat_capacity.Rowlinson_Poling
.. autofunction:: chemicals.heat_capacity.Rowlinson_Bondi
.. autofunction:: chemicals.heat_capacity.Dadgostar_Shaw
.. autofunction:: chemicals.heat_capacity.Dadgostar_Shaw_integral
.. autofunction:: chemicals.heat_capacity.Dadgostar_Shaw_integral_over_T
.. autofunction:: chemicals.heat_capacity.Dadgostar_Shaw_terms

Solid Heat Capacity Estimation Models
-------------------------------------
.. autofunction:: chemicals.heat_capacity.Perry_151
.. autofunction:: chemicals.heat_capacity.Lastovka_solid
.. autofunction:: chemicals.heat_capacity.Lastovka_solid_integral
.. autofunction:: chemicals.heat_capacity.Lastovka_solid_integral_over_T

Utility methods
---------------
.. autoclass:: chemicals.heat_capacity.PiecewiseHeatCapacity

Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: Cp_data_Poling

    Constains data for gases and liquids from [3]_. 
    Simple polynomials for gas heat capacity (not suitable for extrapolation) are available for 308 chemicals. Additionally, constant values in at 298.15 K are available for 348 gases. Constant values in at 298.15 K are available for 245 liquids.

.. data:: TRC_gas_data

    A rigorous expression from [1]_ for modeling gas heat capacity.
    Coefficients for 1961 chemicals are available.

.. data:: CRC_standard_data

    Constant values tabulated in [4]_ at 298.15 K. Data is available for
    533 gases. Data is available for 433 liquids. Data is available for 529
    solids.

.. data:: Cp_dict_PerryI

    Simple polynomials from [5]_ with vaious exponents selected for each expression.
    Coefficients are in units of calories/mol/K. The full expression is
    :math:`C_p = a + bT + c/T^2 + dT^2`. Data is available for 284 compounds.
    Some compounds have gas data, some have liquid data, and have solid 
    (crystal structure) data, sometimes multiple coefficients for different 
    solid phases.

.. data:: zabransky_dicts

    Complicated fits covering different cases and with different forms from [2]_.
    
.. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
    in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
.. [2] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
    Heat Capacity of Liquids: Critical Review and Recommended Values.
    2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
.. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
    New York: McGraw-Hill Professional, 2000.
.. [4] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
    Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
.. [5] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    Eighth Edition. McGraw-Hill Professional, 2007.

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.heat_capacity.Cp_data_Poling

    In [3]: chemicals.heat_capacity.TRC_gas_data

    In [4]: chemicals.heat_capacity.CRC_standard_data

    In [5]: chemicals.heat_capacity.Cp_dict_PerryI['124-38-9'] # gas only

    In [6]: chemicals.heat_capacity.Cp_dict_PerryI['7704-34-9'] # crystal and gas

    In [7]: chemicals.heat_capacity.Cp_dict_PerryI['7440-57-5'] # crystal and liquid

    In [8]: chemicals.heat_capacity.zabransky_dicts.keys()

"""

__all__ = ['heat_capacity_gas_methods',
           'Poling', 'Poling_integral', 'Poling_integral_over_T',
           'Lastovka_Shaw', 'Lastovka_Shaw_integral', 'Lastovka_Shaw_integral_over_T', 
           'Lastovka_Shaw_T_for_Hm', 'Lastovka_Shaw_T_for_Sm', 'Lastovka_Shaw_term_A',
           'TRCCp', 'TRCCp_integral', 'TRCCp_integral_over_T', 
           'heat_capacity_liquid_methods',
           'Rowlinson_Poling', 'Rowlinson_Bondi', 'Dadgostar_Shaw', 
           'Zabransky_quasi_polynomial', 'Zabransky_quasi_polynomial_integral',
           'Zabransky_quasi_polynomial_integral_over_T', 'Zabransky_cubic', 
           'Zabransky_cubic_integral', 'Zabransky_cubic_integral_over_T',
           'Dadgostar_Shaw_integral', 'Dadgostar_Shaw_integral_over_T',
           'Dadgostar_Shaw_terms', 
           'heat_capacity_solid_methods',
           'Lastovka_solid', 'Lastovka_solid_integral', 
           'Lastovka_solid_integral_over_T', 'heat_capacity_solid_methods',
           'ZabranskySpline', 'ZabranskyQuasipolynomial',
           'PiecewiseHeatCapacity',
           ]
import os
from io import open
from chemicals.utils import R, log, exp, to_num, PY37, property_mass_to_molar, source_path, os_path_join, can_load_data
from cmath import log as clog, exp as cexp
from chemicals.data_reader import register_df_source, data_source
from fluids.numerics import newton, brenth, secant, polylog2, numpy as np

__numba_additional_funcs__ = [
    'Lastovka_Shaw_T_for_Hm_err',
    'Lastovka_Shaw_T_for_Sm_err',
]

### Methods introduced in this module

# Gases
TRCIG = 'TRC Thermodynamics of Organic Compounds in the Gas State (1994)'
POLING = 'Poling et al. (2001)'
POLING_CONST = 'Poling et al. (2001) constant'
CRCSTD = 'CRC Standard Thermodynamic Properties of Chemical Substances'
VDI_TABULAR = 'VDI Heat Atlas'
LASTOVKA_SHAW = 'Lastovka and Shaw (2013)'
heat_capacity_gas_methods = (
    TRCIG, POLING, LASTOVKA_SHAW, CRCSTD, POLING_CONST, VDI_TABULAR
)
# Liquids
ZABRANSKY_SPLINE = 'Zabransky spline, averaged heat capacity'
ZABRANSKY_QUASIPOLYNOMIAL = 'Zabransky quasipolynomial, averaged heat capacity'
ZABRANSKY_SPLINE_C = 'Zabransky spline, constant-pressure'
ZABRANSKY_QUASIPOLYNOMIAL_C = 'Zabransky quasipolynomial, constant-pressure'
ZABRANSKY_SPLINE_SAT = 'Zabransky spline, saturation'
ZABRANSKY_QUASIPOLYNOMIAL_SAT = 'Zabransky quasipolynomial, saturation'
ROWLINSON_POLING = 'Rowlinson and Poling (2001)'
ROWLINSON_BONDI = 'Rowlinson and Bondi (1969)'
DADGOSTAR_SHAW = 'Dadgostar and Shaw (2011)'
heat_capacity_liquid_methods = (
    ZABRANSKY_SPLINE, ZABRANSKY_QUASIPOLYNOMIAL, ZABRANSKY_SPLINE_C, 
    ZABRANSKY_QUASIPOLYNOMIAL_C, ZABRANSKY_SPLINE_SAT, ZABRANSKY_QUASIPOLYNOMIAL_SAT,
    VDI_TABULAR, ROWLINSON_POLING, ROWLINSON_BONDI, DADGOSTAR_SHAW, POLING_CONST,
    CRCSTD
)
# Solids
LASTOVKA_S = 'Lastovka, Fulem, Becerra and Shaw (2008)'
PERRY151 = "Perry's Table 2-151"
heat_capacity_solid_methods = (PERRY151, CRCSTD, LASTOVKA_S)

### Heat capacity classes
class ZabranskySpline(object):
    r'''
    Implementation of the cubic spline method presented in [1]_ for 
    calculating the heat capacity of a chemical.
    Implements the enthalpy and entropy integrals as well.
    
    .. math::
        \frac{C}{R}=\sum_{j=0}^3 A_{j+1} \left(\frac{T}{100}\right)^j
    
    Parameters
    ----------
    coeffs : list[float]
        Six coefficients for the equation.
    Tmin : float
        Minimum temperature any experimental data was available at.
    Tmax : float
        Maximum temperature any experimental data was available at.
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    
    '''
    try:
        IS_NUMBA
    except:
        __slots__ = ('coeffs', 'Tmin', 'Tmax')
    
    def __init__(self, coeffs, Tmin, Tmax):
        self.coeffs = coeffs
        self.Tmin = Tmin
        self.Tmax = Tmax

    def calculate(self, T):
        r'''
        Return heat capacity as a function of temperature.
            
        Parameters
        ----------
        T : float
            Temperature, [K]
        
        Returns
        -------
        Cp : float
            Liquid heat capacity as T, [J/mol/K]
        
        '''        
        return Zabransky_cubic(T, *self.coeffs)

    def calculate_integral(self, Ta, Tb):
        r'''
        Return the enthalpy integral of heat capacity from  `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
            
        Returns
        -------
        dS : float
            Enthalpy difference between `Ta` and `Tb`, [J/mol/K]
        
        '''
        return (Zabransky_cubic_integral(Tb, *self.coeffs)
                - Zabransky_cubic_integral(Ta, *self.coeffs))
        
    def calculate_integral_over_T(self, Ta, Tb):
        r'''
        Return the entropy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
            
        Returns
        -------
        dS : float
            Entropy difference between `Ta` and `Tb`, [J/mol/K]
        '''        
        return (Zabransky_cubic_integral_over_T(Tb, *self.coeffs)
                - Zabransky_cubic_integral_over_T(Ta, *self.coeffs))
try:
    if IS_NUMBA:
        ZabranskySpline = jitclass([('coeffs', numba.types.UniTuple(numba.float64, 4)),
                ('Tmin', numba.float64),
                ('Tmax', numba.float64)])(ZabranskySpline)
except:
    pass

class ZabranskyQuasipolynomial(object):
    r'''
    Quasi-polynomial object for calculating the heat capacity of a chemical.
    Implements the enthalpy and entropy integrals as well.
    
    .. math::
        \frac{C}{R}=A_1\ln(1-T_r) + \frac{A_2}{1-T_r}
        + \sum_{j=0}^m A_{j+3} T_r^j
    
    Parameters
    ----------
    coeffs : list[float]
        Six coefficients for the equation.
    Tc : float
        Critical temperature of the chemical, as used in the formula.
    Tmin : float
        Minimum temperature any experimental data was available at.
    Tmax : float
        Maximum temperature any experimental data was available at.
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    
    '''
    try:
        IS_NUMBA
    except:
        __slots__ = ('coeffs', 'Tc', 'Tmin', 'Tmax')
    
    def __init__(self, coeffs, Tc, Tmin, Tmax):
        self.coeffs = coeffs
        self.Tc = Tc
        self.Tmin = Tmin
        self.Tmax = Tmax

    def calculate(self, T):
        r'''
        Return the heat capacity as a function of temperature.
            
        Parameters
        ----------
        T : float
            Temperature, [K]
        
        Returns
        -------
        Cp : float
            Liquid heat capacity as T, [J/mol/K]
        
        '''
        return Zabransky_quasi_polynomial(T, self.Tc, *self.coeffs)
                                          
    def calculate_integral(self, Ta, Tb):
        r'''
        Return the enthalpy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
            
        Returns
        -------
        dH : float
            Enthalpy difference between `Ta` and `Tb`, [J/mol]
        
        '''        
        return (Zabransky_quasi_polynomial_integral(Tb, self.Tc, *self.coeffs)
               - Zabransky_quasi_polynomial_integral(Ta, self.Tc, *self.coeffs))
    
    def calculate_integral_over_T(self, Ta, Tb):
        r'''
        Return the entropy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
        
        Returns
        -------
        dS : float
            Entropy difference between `Ta` and `Tb`, [J/mol/K]
        
        '''        
        return (Zabransky_quasi_polynomial_integral_over_T(Tb, self.Tc, *self.coeffs)
               - Zabransky_quasi_polynomial_integral_over_T(Ta, self.Tc, *self.coeffs))
try:
    if IS_NUMBA:
        ZabranskyQuasipolynomial = jitclass([('coeffs', numba.types.UniTuple(numba.float64, 6)),
            ('Tc', numba.float64),
            ('Tmin', numba.float64),
            ('Tmax', numba.float64)])(ZabranskyQuasipolynomial)
except:
    pass


class PiecewiseHeatCapacity(object):
    r"""
    Create a PiecewiseHeatCapacity object for calculating heat capacity and the 
    enthalpy and entropy integrals using piecewise models.
    
    Parameters
    ----------
    models : Iterable[HeatCapacity]
        Piecewise heat capacity objects.
    
    """
    # Dev note - not possible to jitclass this as the model types are not explicit
    __slots__ = ('models', 'Tmin', 'Tmax')
    
    def __init__(self, models):
        self.models = tuple(sorted(models, key=lambda x: x.Tmin))
        self.Tmin = self.models[0].Tmin
        self.Tmax = self.models[-1].Tmax
    
    def __iter__(self):
        return self.models.__iter__()
    
    def calculate(self, T):
        r'''
        Return the heat capacity as a function of temperature.
            
        Parameters
        ----------
        T : float
            Temperature, [K]
        
        Raises
        ------
        ValueError
            If the temperature in not within the domain of any of the models 
            (i.e if Tmin <= T <= Tmax cannot be satisfied by any of the models).
        
        See also
        --------
        PiecewiseHeatCapacity.force_calculate
        
        Returns
        -------
        Cp : float
            Liquid heat capacity as T, [J/mol/K]
        
        '''
        if T >= self.Tmin:
            for model in self.models:
                if T <= model.Tmax: return model.calculate(T)
        raise ValueError("no valid model at T=%g K" % T)
    
    def force_calculate(self, T):
        r'''
        Return the heat capacity as a function of temperature. 
            
        Parameters
        ----------
        T : float
            Temperature, [K]
        
        Notes
        -----
        This method extrapolates when temperature is not within the domain of 
        any of the models (i.e if Tmin <= T <= Tmax cannot be satisfied by any
        of the models).
        
        See also
        --------
        PiecewiseHeatCapacity.calculate
        
        Returns
        -------
        Cp : float
            Liquid heat capacity as T, [J/mol/K]
        
        '''
        for model in self.models:
            if T <= model.Tmax: break
        return model.calculate(T)
    
    def calculate_integral(self, Ta, Tb):
        r'''
        Return the enthalpy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
        
        Raises
        ------
        ValueError
            If the temperature in not within the domain of any of the models 
            (i.e if Tmin <= T <= Tmax cannot be satisfied by any of the models).
        
        Notes
        -----
        Analytically integrates piecewise through all models.
        
        See also
        --------
        PiecewiseHeatCapacity.force_calculate_integral
        
        Returns
        -------
        dH : float
            Enthalpy difference between `Ta` and `Tb`, [J/mol]
        
        '''
        if Tb < Ta: return -self.calculate_integral(Tb, Ta)
        if Ta < self.Tmin: raise ValueError("no valid model at T=%g K" % Ta)
        elif Tb > self.Tmax: raise ValueError("no valid model at T=%g K" % Tb)
        return self.force_calculate_integral(Ta, Tb)
    
    def force_calculate_integral(self, Ta, Tb):
        r'''
        Return the enthalpy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
        
        Raises
        ------
        ValueError
            If the temperature in not within the domain of any of the models 
            (i.e if Tmin <= T <= Tmax cannot be satisfied by any of the models).
        
        Notes
        -----
        Analytically integrates piecewise through all models and extrapolates
        when temperature is not within the domain of any of the models 
        (i.e if Tmin <= T <= Tmax cannot be satisfied by any of the models).
        
        See also
        --------
        PiecewiseHeatCapacity.calculate_integral
        
        Returns
        -------
        dH : float
            Enthalpy difference between `Ta` and `Tb`, [J/mol]
        
        '''
        if Tb < Ta: return -self.force_calculate_integral(Tb, Ta)
        integral = 0.
        for model in self.models:
            Tmax = model.Tmax
            if Tb <= Tmax:
                return integral + model.calculate_integral(Ta, Tb)
            else:
                integral += model.calculate_integral(Ta, Tmax)
                Ta = Tmax
        return integral + model.calculate_integral(Ta, Tb)
    
    def calculate_integral_over_T(self, Ta, Tb):
        r'''
        Return the entropy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
            
        Notes
        -----
        Analytically integrates piecewise through all models.
        
        Raises
        ------
        ValueError
            If the temperature in not within the domain of any of the models 
            (i.e if Tmin <= T <= Tmax cannot be satisfied by any of the models).
        
        See also
        --------
        PiecewiseHeatCapacity.force_calculate_integral_over_T
        
        Returns
        -------
        dS : float
            Entropy difference between `Ta` and `Tb`, [J/mol/K]
        
        '''
        if Tb < Ta: return -self.calculate_integral_over_T(Tb, Ta)
        if Ta < self.Tmin: raise ValueError("no valid model at T=%d K" % Ta)
        elif Tb > self.Tmax: raise ValueError("no valid model at T=%d K" % Tb)
        return self.force_calculate_integral_over_T(Ta, Tb)
        
    def force_calculate_integral_over_T(self, Ta, Tb):
        r'''
        Return the entropy integral of heat capacity from `Ta` to `Tb`.
            
        Parameters
        ----------
        Ta : float
            Initial temperature, [K]
        Tb : float
            Final temperature, [K]
            
        Notes
        -----
        Analytically integrates piecewise through all models and extrapolates
        when temperature is not within the domain of any of the models 
        (i.e if Tmin <= T <= Tmax cannot be satisfied by any of the models).
        
        See also
        --------
        PiecewiseHeatCapacity.calculate_integral_over_T
        
        Returns
        -------
        dS : float
            Entropy difference between `Ta` and `Tb`, [J/mol/K]
        
        '''
        if Tb < Ta: return -self.force_calculate_integral_over_T(Tb, Ta)
        integral = 0.
        for model in self.models:
            Tmax = model.Tmax
            if Tb <= Tmax:
                return integral + model.calculate_integral_over_T(Ta, Tb)
            else:
                integral += model.calculate_integral_over_T(Ta, Tmax)
                Ta = Tmax
        return integral + model.calculate_integral_over_T(Ta, Tb)


### Register data sources and lazy load them

folder = os_path_join(source_path, 'Heat Capacity')
register_df_source(folder, 'PolingDatabank.tsv')
register_df_source(folder, 'TRC Thermodynamics of Organic Compounds in the Gas State.tsv', csv_kwargs={
    'dtype':{'Tmin': float, 'Tmax': float, 'a0': float, 'a1': float, 'a2': float,
             'a3': float, 'a4': float, 'a5': float, 'a6': float, 'a7': float, 
             'I': float, 'J': float, 'Hfg': float}})

register_df_source(folder, 'CRC Standard Thermodynamic Properties of Chemical Substances.tsv')

_Cp_data_loaded = False
def _load_Cp_data():
    global Cp_data_Poling, Cp_values_Poling, TRC_gas_data, TRC_gas_values
    global CRC_standard_data, Cp_dict_PerryI
    global zabransky_dict_sat_s, zabransky_dict_sat_p, zabransky_dict_const_s
    global zabransky_dict_const_p, zabransky_dict_iso_s, zabransky_dict_iso_p
    global type_to_zabransky_dict, zabransky_dicts, _Cp_data_loaded
    Cp_data_Poling = data_source('PolingDatabank.tsv')
    TRC_gas_data = data_source('TRC Thermodynamics of Organic Compounds in the Gas State.tsv')
    CRC_standard_data = data_source('CRC Standard Thermodynamic Properties of Chemical Substances.tsv')
    
    TRC_gas_values = np.array(TRC_gas_data.values[:, 1:], dtype=float)
    Cp_values_Poling = np.array(Cp_data_Poling.values[:, 1:], dtype=float)
    
    # Read in a dict of heat capacities of irnorganic and elemental solids.
    # These are in section 2, table 151 in:
    # Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    # Eighth Edition. McGraw-Hill Professional, 2007.
    # Formula:
    # Cp(Cal/mol/K) = Const + Lin*T + Quadinv/T^2 + Quadinv*T^2
    # Phases:
    # c, gls, l, g.
    zabransky_dict_sat_s = {}
    zabransky_dict_sat_p = {}
    zabransky_dict_const_s = {}
    zabransky_dict_const_p = {}
    zabransky_dict_iso_s = {}
    zabransky_dict_iso_p = {}
    # C means average heat capacity values, from less rigorous experiments
    # sat means heat capacity along the saturation line
    # p means constant-pressure values, 
    # second argument is whether or not it has a spline
    type_to_zabransky_dict = {
        ('C', True): zabransky_dict_const_s, 
        ('C', False):   zabransky_dict_const_p,
        ('sat', True):  zabransky_dict_sat_s,
        ('sat', False): zabransky_dict_sat_p,
        ('p', True):    zabransky_dict_iso_s,
        ('p', False):   zabransky_dict_iso_p
    }
    zabransky_dicts = {
        ZABRANSKY_SPLINE: zabransky_dict_const_s,
        ZABRANSKY_QUASIPOLYNOMIAL: zabransky_dict_const_p,
        ZABRANSKY_SPLINE_C: zabransky_dict_iso_s,
        ZABRANSKY_QUASIPOLYNOMIAL_C: zabransky_dict_iso_p,
        ZABRANSKY_SPLINE_SAT: zabransky_dict_sat_s,
        ZABRANSKY_QUASIPOLYNOMIAL_SAT: zabransky_dict_sat_p
    }
    with open(os.path.join(folder, 'Zabransky.tsv'), encoding='utf-8') as f:
        next(f)
        for line in f:
            values = to_num(line.strip('\n').split('\t'))
            (CAS, name, Type, uncertainty, Tmin, Tmax,
             a1s, a2s, a3s, a4s, a1p, a2p, a3p, a4p, a5p, a6p, Tc) = values
            spline = bool(a1s) # False if Quasypolynomial, True if spline
            d = type_to_zabransky_dict[(Type, spline)]
            if spline:
                coeffs = (a1s, a2s, a3s, a4s)
                if CAS not in d:
                    d[CAS] = [ZabranskySpline(coeffs, Tmin, Tmax)]
                else:
                    d[CAS].append(ZabranskySpline(coeffs, Tmin, Tmax))
            else:
                # No duplicates for quasipolynomials
                coeffs = (a1p, a2p, a3p, a4p, a5p, a6p)
                d[CAS] = ZabranskyQuasipolynomial(coeffs, Tc, Tmin, Tmax)
    for dct in (zabransky_dict_const_s, zabransky_dict_iso_s, zabransky_dict_sat_s):
        for CAS in dct: dct[CAS] = PiecewiseHeatCapacity(dct[CAS])
    # Used to generate data. Do not delete!
    # Cp_dict_PerryI = {}
    # with open(os.path.join(folder, 'Perrys Table 2-151.tsv'), encoding='utf-8') as f:
    #     '''Read in a dict of heat capacities of irnorganic and elemental solids.
    #     These are in section 2, table 151 in:
    #     Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    #     Eighth Edition. McGraw-Hill Professional, 2007.
    #     Formula:
    #     Cp(Cal/mol/K) = Const + Lin*T + Quadinv/T^2 + Quadinv*T^2
    #     Phases: c, gls, l, g.
    #     '''
    #     next(f)
    #     for line in f:
    #         values = to_num(line.strip('\n').split('\t'))
    #         (CASRN, _formula, _phase, _subphase, Const, Lin, Quadinv, Quad, Tmin,
    #          Tmax, err) = values
    #         if Lin is None:
    #             Lin = 0
    #         if Quadinv is None:
    #             Quadinv = 0
    #         if Quad is None:
    #             Quad = 0
    #         if CASRN in _PerryI and CASRN:
    #             a = _PerryI[CASRN]
    #             a.update({_phase: {"Formula": _formula, "Phase": _phase,
    #                                "Subphase": _subphase, "Const": Const,
    #                                "Lin": Lin, "Quadinv": Quadinv, "Quad": Quad,
    #                                "Tmin": Tmin, "Tmax": Tmax, "Error": err}})
    #             _PerryI[CASRN] = a
    #         else:
    #             _PerryI[CASRN] = {_phase: {"Formula": _formula, "Phase": _phase,
    #                                        "Subphase": _subphase, "Const": Const,
    #                                        "Lin": Lin, "Quadinv": Quadinv,
    #                                        "Quad": Quad, "Tmin": Tmin,
    #                                        "Tmax": Tmax, "Error": err}}
    '''
    Read in a dict of 2481 thermodynamic property sets of different phases from:
    Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
    Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    Warning: 11 duplicated chemicals are present and currently clobbered.
    '''
    import json
    with open(os.path.join(folder, 'Perrys Table 2-151.json')) as f:
        Cp_dict_PerryI = json.loads(f.read())
    _Cp_data_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('Cp_data_Poling', 'Cp_values_Poling', 'TRC_gas_data', 'TRC_gas_values', 'CRC_standard_data',
                    'Cp_dict_PerryI', 'zabransky_dict_sat_s', 'zabransky_dict_sat_p', 
                    'zabransky_dict_const_s', 'zabransky_dict_const_p', 'zabransky_dict_iso_s',  
                    'zabransky_dict_iso_p', 'type_to_zabransky_dict', 'zabransky_dicts'):
            _load_Cp_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_Cp_data()



### Heat capacities of gases

def Poling(T, a, b, c, d, e):
    r"""
    Return the ideal-gas molar heat capacity of a chemical using 
    polynomial regressed coefficients as described by Poling et. al. [1]_.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    a,b,c,d,e : float
        Regressed coefficients.
    
    Returns
    -------
    Cpgm : float
        Gas molar heat capacity, [J/mol/K]
    
    Notes
    -----
    The ideal gas heat capacity is given by:
    
    .. math:: C_n = R*(a + bT + cT^2 + dT^3 + eT^4)
    
    The data is based on the Poling data bank.
    
    See Also
    --------
    Poling_integral
    Poling_integral_over_T
    
    Examples
    --------
    Compute the gas heat capacity of Methane at 300 K:
    
    >>> Poling(T=300., a=4.568, b=-0.008975, c=3.631e-05, d=-3.407e-08, e=1.091e-11)
    35.850973388425
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
           New York: McGraw-Hill Professional, 2000.
    
    """
    return R*(T*(T*(T*(T*e + d) + c) + b) + a)
 
def Poling_integral(T, a, b, c, d, e):
    r"""
    Return the integral of the ideal-gas constant-pressure heat capacity 
    of a chemical using polynomial regressed coefficients as described by
    Poling et. al. [1]_.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    a,b,c,d,e : float
        Regressed coefficients.
    
    Returns
    -------
    H : float
        Difference in enthalpy from 0 K, [J/mol]
        
    Notes
    -----
    Integral was computed with SymPy.
    
    See Also
    --------
    Poling
    Poling_integral_over_T
    
    Examples
    --------
    Compute the gas enthalpy of Methane at 300 K (with reference to 0 K):
    
    >>> Poling_integral(T=300., a=4.568, b=-0.008975, c=3.631e-05, d=-3.407e-08, e=1.091e-11)
    10223.67533722261
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
           New York: McGraw-Hill Professional, 2000.
    
    """
    return R*(((((0.2*e)*T + 0.25*d)*T + c*(1.0/3.))*T + 0.5*b)*T + a)*T
 
def Poling_integral_over_T(T, a, b, c, d, e):
    r"""
    Return the integral over temperature of the ideal-gas constant-pressure 
    heat capacity of a chemical using polynomial regressed coefficients as
    described by Poling et. al. [1]_.
    
    Parameters
    ----------
    T : float
        Temperature, [K]
    a,b,c,d,e : float
        Regressed coefficients.
    
    Returns
    -------
    S : float
        Difference in entropy from 0 K, [J/mol/K]
        
    Notes
    -----
    Integral was computed with SymPy.
    
    See Also
    --------
    Poling
    Poling_integral
    
    Examples
    --------
    Compute the gas entropy of Methane at 300 K (with reference to 0 K):
    
    >>> Poling_integral_over_T(T=300., a=4.568, b=-0.008975, c=3.631e-05, d=-3.407e-08, e=1.091e-11)
    205.46526328058
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
           New York: McGraw-Hill Professional, 2000.
    
    """
    return R*(((((0.25*e)*T + d*(1.0/3.))*T + 0.5*c)*T + b)*T + a*log(T))

def Lastovka_Shaw_term_A(similarity_variable, cyclic_aliphatic):
    """
    Return Term A in Lastovka-Shaw equation.

    Parameters
    ----------
    similarity_variable : float
        Similarity variable as defined in [1]_, [mol/g]
    cyclic_aliphatic: bool, optional
        Whether or not chemical is cyclic aliphatic, [-]

    Returns
    -------
    term_A : float
        Term A in Lastovka-Shaw equation, [J/g]

    See Also
    --------
    Lastovka_Shaw
    Lastovka_Shaw_integral
    Lastovka_Shaw_integral_over_T

    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.

    """
    a = similarity_variable
    if cyclic_aliphatic:
        A1 = -0.1793547
        A2 = 3.86944439
        term_A = A1 + A2*a
    else:
        # A1 = 0.58
        A2 = 1.25
        A1_minus_A2 = -0.67 # (A1 - A2)
        A3 = 0.17338003 # 803 instead of 8003 in another paper
        # A4 = 0.014
        A4_inv = 71.42857142857143 # 1 / A4
        term_A = A2 + A1_minus_A2/(1. + exp((a-A3)*A4_inv)) # One reference says exp((a-A3)/A4)
        # Personal communication confirms the change
    return term_A

def Lastovka_Shaw(T, similarity_variable, cyclic_aliphatic=False, MW=None, term_A=None):
    r'''Calculate ideal-gas constant-pressure heat capacity with the similarity
    variable concept and method as shown in [1]_.
    
    .. math::
        term_A = A1 + A2*a \text{ if cyclic aliphatic} 
        
    .. math::
        term_A = \left(A_2 + \frac{A_1 - A_2}{1 + \exp(\frac{\alpha-A_3}{A_4})}\right) \text{ if not cyclic aliphatic}
        
    .. math::
        C_p^0 = term_A
        + (B_{11} + B_{12}\alpha)\left(-\frac{(C_{11} + C_{12}\alpha)}{T}\right)^2
        \frac{\exp(-(C_{11} + C_{12}\alpha)/T)}{[1-\exp(-(C_{11}+C_{12}\alpha)/T)]^2}
        + (B_{21} + B_{22}\alpha)\left(-\frac{(C_{21} + C_{22}\alpha)}{T}\right)^2
        \frac{\exp(-(C_{21} + C_{22}\alpha)/T)}{[1-\exp(-(C_{21}+C_{22}\alpha)/T)]^2}
        
    Parameters
    ----------
    T : float
        Temperature of gas [K]
    similarity_variable : float
        Similarity variable as defined in [1]_, [mol/g]
    cyclic_aliphatic: bool, optional
        Whether or not chemical is cyclic aliphatic, [-]
    MW : float, optional
        Molecular weight, [g/mol]
    term_A : float, optional
        Term A in Lastovka-Shaw equation, [J/g]
    
    Returns
    -------
    Cpg : float
        Gas constant-pressure heat capacity, [J/mol/K if MW given; J/kg/K otherwise]
        
    Notes
    -----
    Original model is in terms of J/g/K.
    
    A1 = -0.1793547 \text{ if cyclic aliphatic}
    
    A1 = 0.58 \text{ if not cyclic aliphatic}
    
    A2 = 3.86944439 \text{ if cyclic aliphatic}
    
    A2 = 1.25 \text{ if not cyclic aliphatic}
    
    A3 = 0.17338003
    
    A4 = 0.014
    
    B11 = 0.73917383
    
    B12 = 8.88308889
    
    C11 = 1188.28051
    
    C12 = 1813.04613
    
    B21 = 0.0483019
    
    B22 = 4.35656721
    
    C21 = 2897.01927
    
    C22 = 5987.80407
    
    Examples
    --------
    Estimate the heat capacity of n-decane gas in J/kg/K:
    
    >>> Lastovka_Shaw(1000.0, 0.22491)
    3730.2807601773725
    
    Estimate the heat capacity of n-decane gas in J/mol/K:
    
    >>> Lastovka_Shaw(1000.0, 0.22491, MW=142.28)
    530.7443465580366
    
    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    
    '''
    a = similarity_variable
    if term_A is None: term_A = Lastovka_Shaw_term_A(a, cyclic_aliphatic)

    T_inv = 1.0/T
    B11 = 0.73917383
    B12 = 8.88308889
    C11 = 1188.28051
    C12 = 1813.04613
    B21 = 0.0483019
    B22 = 4.35656721
    C21 = 2897.01927
    C22 = 5987.80407
    C11_C12a_T =(C11+C12*a)*T_inv
    expm_C11_C12a_T = exp(-C11_C12a_T)
    x1 = 1.0/(1.0 - expm_C11_C12a_T)
    
    C21_C22a_T = (C21+C22*a)*T_inv
    expm_C21_C22a_T = exp(-C21_C22a_T)
    x2 = 1.0/(1.0 - expm_C21_C22a_T)
    
    Cp = term_A + (B11 + B12*a)*(C11_C12a_T*C11_C12a_T)*expm_C11_C12a_T*x1*x1
    Cp += (B21 + B22*a)*(C21_C22a_T*C21_C22a_T)*expm_C21_C22a_T*x2*x2
    return Cp*1000. if MW is None else Cp*MW

def Lastovka_Shaw_integral(T, similarity_variable, cyclic_aliphatic=False,
                           MW=None, term_A=None):
    r'''Calculate the integral of ideal-gas constant-pressure heat capacity 
    with the similarity variable concept and method as shown in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature of gas [K]
    cyclic_aliphatic: bool, optional
        Whether or not chemical is cyclic aliphatic, [-]
    MW : float, optional
        Molecular weight, [g/mol]
    term_A : float, optional
        Term A in Lastovka-Shaw equation, [J/g]
        
    Returns
    -------
    H : float
        Difference in enthalpy from 0 K, [J/mol if MW given; J/kg otherwise]
        
    Notes
    -----
    Original model is in terms of J/g/K. 
    Integral was computed with SymPy.
    
    See Also
    --------
    Lastovka_Shaw
    Lastovka_Shaw_integral_over_T
    
    Examples
    --------
    >>> Lastovka_Shaw_integral(300.0, 0.1333)
    5283095.816018478
    
    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    
    '''
    a = similarity_variable
    if term_A is None: term_A = Lastovka_Shaw_term_A(a, cyclic_aliphatic)

    B11 = 0.73917383
    B12 = 8.88308889
    C11 = 1188.28051
    C12 = 1813.04613
    B21 = 0.0483019
    B22 = 4.35656721
    C21 = 2897.01927
    C22 = 5987.80407
    x1 = -C11 - C12*a
    x2 = -C21 - C22*a
    T_inv = 1.0/T
    
    H = (T*term_A - (B11 + B12*a)*(x1*x1)/(x1 - x1*exp(x1*T_inv)) 
                  - (B21 + B22*a)*(x2*x2)/(x2 - x2*exp(x2*T_inv)))
    return H*1000. if MW is None else H*MW

def Lastovka_Shaw_integral_over_T(T, similarity_variable, cyclic_aliphatic=False,
                                  MW=None, term_A=None):
    r'''Calculate the integral over temperature of ideal-gas constant-pressure 
    heat capacity with the similarity variable concept and method as shown in
    [1]_.
    
    Parameters
    ----------
    T : float
        Temperature of gas [K]
    similarity_variable : float
        Similarity variable as defined in [1]_, [mol/g]
    cyclic_aliphatic: bool, optional
        Whether or not chemical is cyclic aliphatic, [-]
    MW : float, optional
        Molecular weight, [g/mol]
    term_A : float, optional
        Term A in Lastovka-Shaw equation, [J/g]
        
    Returns
    -------
    S : float
        Difference in entropy from 0 K, [J/mol/K if MW given; J/kg/K otherwise]
        
    Notes
    -----
    Original model is in terms of J/g/K. Note that the model is for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    Integral was computed with SymPy.
    
    See Also
    --------
    Lastovka_Shaw
    Lastovka_Shaw_integral
    
    Examples
    --------
    >>> Lastovka_Shaw_integral_over_T(300.0, 0.1333)
    3609.791928945323
    
    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    '''
    a = similarity_variable
    if term_A is None: term_A = Lastovka_Shaw_term_A(a, cyclic_aliphatic)

    T_inv = 1.0/T
    a2 = a*a
    B11 = 0.73917383
    B12 = 8.88308889
    C11 = 1188.28051
    C12 = 1813.04613
    B21 = 0.0483019
    B22 = 4.35656721
    C21 = 2897.01927
    C22 = 5987.80407
    S = (term_A*clog(T) + (-B11 - B12*a)*clog(cexp((-C11 - C12*a)*T_inv) - 1.) 
        + (-B11*C11 - B11*C12*a - B12*C11*a - B12*C12*a2)/(T*cexp((-C11
        - C12*a)*T_inv) - T) - (B11*C11 + B11*C12*a + B12*C11*a + B12*C12*a2)*T_inv)
    S += ((-B21 - B22*a)*clog(cexp((-C21 - C22*a)*T_inv) - 1.) + (-B21*C21 - B21*C22*a
        - B22*C21*a - B22*C22*a2)/(T*cexp((-C21 - C22*a)*T_inv) - T) - (B21*C21
        + B21*C22*a + B22*C21*a + B22*C22*a2)*T_inv)
    # There is a non-real component, but it is only a function of similariy 
    # variable and so will always cancel out.
    return S.real*1000. if MW is None else S.real*MW

def Lastovka_Shaw_T_for_Hm_err(T, MW, similarity_variable, H_ref, Hm, cyclic_aliphatic, term_A):
    H1 = Lastovka_Shaw_integral(T, similarity_variable, cyclic_aliphatic, MW, term_A)
    dH = H1 - H_ref
    err = (dH - Hm)
    return err

def Lastovka_Shaw_T_for_Hm(Hm, MW, similarity_variable, T_ref=298.15, 
                           factor=1.0, cyclic_aliphatic=None, term_A=None):
    r'''Uses the Lastovka-Shaw ideal-gas heat capacity correlation to solve for
    the temperature which has a specified `Hm`, as is required in PH flashes,
    as shown in [1]_.

    Parameters
    ----------
    Hm : float
        Molar enthalpy spec, [J/mol]
    MW : float
        Molecular weight of the pure compound or mixture average, [g/mol]
    similarity_variable : float
        Similarity variable as defined in [1]_, [mol/g]
    T_ref : float, optional
        Reference enthlapy temperature, [K]
    factor : float, optional
        A factor to increase or decrease the predicted value of the 
        method, [-]
    cyclic_aliphatic: bool, optional
        Whether or not chemical is cyclic aliphatic, [-]
    term_A : float, optional
        Term A in Lastovka-Shaw equation, [J/g]

    Returns
    -------
    T : float
        Temperature of gas to meet the molar enthalpy spec, [K]

    Notes
    -----

    See Also
    --------
    Lastovka_Shaw
    Lastovka_Shaw_integral
    Lastovka_Shaw_integral_over_T

    Examples
    --------
    >>> Lastovka_Shaw_T_for_Hm(Hm=55000, MW=80.0, similarity_variable=0.23)
    600.0943429567602
    
    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    '''
    Hm /= factor
    a = similarity_variable
    if term_A is None: term_A = Lastovka_Shaw_term_A(a, cyclic_aliphatic)
    H_ref = Lastovka_Shaw_integral(T_ref, similarity_variable, cyclic_aliphatic, MW, term_A)
    args = (MW, a, H_ref, Hm, cyclic_aliphatic, term_A)
    try:
        return secant(Lastovka_Shaw_T_for_Hm_err, 500.0, ytol=1e-4, args=args)
    except:
        try:
            return brenth(Lastovka_Shaw_T_for_Hm_err, 1e-3, 1e5, args)
        except:
            if Lastovka_Shaw_T_for_Hm_err(1e-11, *args) > 0:
                raise ValueError("For gas only enthalpy spec to be correct, "
                                 "model requires negative temperature")
            raise ValueError("Could not converge")

def Lastovka_Shaw_T_for_Sm_err(T, MW, similarity_variable, S_ref, Sm, cyclic_aliphatic, term_A):
    S1 = Lastovka_Shaw_integral_over_T(T, similarity_variable, cyclic_aliphatic, MW, term_A)
    dS = S1 - S_ref
    err = (dS - Sm)
    return err

def Lastovka_Shaw_T_for_Sm(Sm, MW, similarity_variable, T_ref=298.15, 
                           factor=1.0, cyclic_aliphatic=None, term_A=None):
    r'''Uses the Lastovka-Shaw ideal-gas heat capacity correlation to solve for
    the temperature which has a specified `Sm`, as is required in PS flashes,
    as shown in [1]_.

    Parameters
    ----------
    Sm : float
        Molar entropy spec, [J/mol/K]
    MW : float
        Molecular weight of the pure compound or mixture average, [g/mol]
    similarity_variable : float
        Similarity variable as defined in [1]_, [mol/g]
    T_ref : float, optional
        Reference enthlapy temperature, [K]
    factor : float, optional
        A factor to increase or decrease the predicted value of the 
        method, [-]
    cyclic_aliphatic: bool, optional
        Whether or not chemical is cyclic aliphatic, [-]
    term_A : float, optional
        Term A in Lastovka-Shaw equation, [J/g]

    Returns
    -------
    T : float
        Temperature of gas to meet the molar entropy spec, [K]

    Notes
    -----

    See Also
    --------
    Lastovka_Shaw
    Lastovka_Shaw_integral
    Lastovka_Shaw_integral_over_T

    Examples
    --------
    >>> Lastovka_Shaw_T_for_Sm(Sm=112.80, MW=72.151, similarity_variable=0.2356)
    603.4298291570276
    
    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    '''
    Sm /= factor
    a = similarity_variable
    if term_A is None: term_A = Lastovka_Shaw_term_A(a, cyclic_aliphatic)
    S_ref = Lastovka_Shaw_integral_over_T(T_ref, a, cyclic_aliphatic, MW, term_A)
    args = (MW, a, S_ref, Sm, cyclic_aliphatic, term_A)
    try:
        return secant(Lastovka_Shaw_T_for_Sm_err, 500, ytol=1e-4, high=10000,
                      args=args)
    except:
        try:
            return brenth(Lastovka_Shaw_T_for_Sm_err, 1e-3, 1e5,
                          args=args)
        except:
            if Lastovka_Shaw_T_for_Sm_err(1e-11, *args) > 0.0:
                raise ValueError("For gas only entropy spec to be correct, "
                                 "model requires negative temperature")
            raise ValueError("Could not converge")

def TRCCp(T, a0, a1, a2, a3, a4, a5, a6, a7):
    r'''Calculates ideal gas heat capacity using the model developed in [1]_.
    The ideal gas heat capacity is given by:
        
    .. math::
        C_p = R\left(a_0 + (a_1/T^2) \exp(-a_2/T) + a_3 y^2
        + (a_4 - a_5/(T-a_7)^2 )y^j \right)
        
    .. math::
        y = \frac{T-a_7}{T+a_6} \text{ for } T > a_7 \text{ otherwise } 0
        
    Parameters
    ----------
    T : float
        Temperature [K]
    a1-a7 : float
        Coefficients
        
    Returns
    -------
    Cp : float
        Ideal gas heat capacity , [J/mol/K]
        
    Notes
    -----
    j is set to 8. Analytical integrals are available for this expression.
    
    Examples
    --------
    >>> TRCCp(300, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201.)
    42.065271080974654
    
    References
    ----------
    .. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
       in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
    '''
    if T <= a7:
        y = 0.
    else:
        y = (T - a7)/(T + a6)
    T_inv = 1.0/T
    y2 = y*y
    T_m_a7 = T - a7
    if T == a7:
        # When T_m_a7 approaches 0, T = a7
        Cp = R*(a0 + (a1*T_inv*T_inv)*exp(-a2*T_inv) + y2*(a3 + (a4)*y2*y2*y2))
    else:
        Cp = R*(a0 + (a1*T_inv*T_inv)*exp(-a2*T_inv) + y2*(a3 + (a4 - a5/(T_m_a7*T_m_a7))*y2*y2*y2))
    return Cp

def TRCCp_integral(T, a0, a1, a2, a3, a4, a5, a6, a7, I=0):
    r'''Integrates ideal gas heat capacity using the model developed in [1]_.
    Best used as a delta only.
    The difference in enthalpy with respect to 0 K is given by:
        
    .. math::
        \frac{H(T) - H^{ref}}{RT} = a_0 + a_1x(a_2)/(a_2T) + I/T + h(T)/T
        
    .. math::
        h(T) = (a_5 + a_7)\left[(2a_3 + 8a_4)\ln(1-y)+ \left\{a_3\left(1 + 
        \frac{1}{1-y}\right) + a_4\left(7 + \frac{1}{1-y}\right)\right\}y
        + a_4\left\{3y^2 + (5/3)y^3 + y^4 + (3/5)y^5 + (1/3)y^6\right\} 
        + (1/7)\left\{a_4 - \frac{a_5}{(a_6+a_7)^2}\right\}y^7\right]
        
    .. math::
        h(T) = 0 \text{ for } T \le a_7
        y = \frac{T-a_7}{T+a_6} \text{ for } T > a_7 \text{ otherwise } 0
        
    Parameters
    ----------
    T : float
        Temperature [K]
    a1-a7 : float
        Coefficients
    I : float, optional
        Integral offset
        
    Returns
    -------
    H-H(0) : float
        Difference in enthalpy from 0 K , [J/mol]
        
    Notes
    -----
    Analytical integral as provided in [1]_ and verified with numerical
    integration. 
    
    Examples
    --------
    >>> TRCCp_integral(298.15, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 
    ... 201., 1.2)
    10802.536262068483
    
    References
    ----------
    .. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
       in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
    '''
    if T <= a7:
        y = 0.
    else:
        y = (T - a7)/(T + a6)
    y2 = y*y
    y4 = y2*y2
    if T <= a7:
        h = 0.0
    else:
        first = a6 + a7
        one_m_y = 1.0 - y
        second = (2.*a3 + 8.*a4)*log(one_m_y)
        third = (a3*(1. + 1./(one_m_y)) + a4*(7. + 1./(one_m_y)))*y
        fourth = a4*(3.*y2 + 5./3.*y*y2 + y4 + 0.6*y4*y + 1/3.*y4*y2)
        fifth = 1/7.*(a4 - a5/(first*first))*y4*y2*y
        h = first*(second + third + fourth + fifth)
    return (a0 + a1*exp(-a2/T)/(a2*T) + I/T + h/T)*R*T

def TRCCp_integral_over_T(T, a0, a1, a2, a3, a4, a5, a6, a7, J=0):
    r'''Integrates ideal gas heat capacity over T using the model developed in 
    [1]_. Best used as a delta only.
    The difference in ideal-gas entropy with respect to 0 K is given by:
        
    .. math::
        \frac{S^\circ}{R} = J + a_0\ln T + \frac{a_1}{a_2^2}\left(1
        + \frac{a_2}{T}\right)x(a_2) + s(T)
        s(T) = \left[\left\{a_3 + \left(\frac{a_4 a_7^2 - a_5}{a_6^2}\right)
        \left(\frac{a_7}{a_6}\right)^4\right\}\left(\frac{a_7}{a_6}\right)^2
        \ln z + (a_3 + a_4)\ln\left(\frac{T+a_6}{a_6+a_7}\right)
        +\sum_{i=1}^7 \left\{\left(\frac{a_4 a_7^2 - a_5}{a_6^2}\right)\left(
        \frac{-a_7}{a_6}\right)^{6-i} - a_4\right\}\frac{y^i}{i}
        - \left\{\frac{a_3}{a_6}(a_6 + a_7) + \frac{a_5 y^6}{7a_7(a_6+a_7)}
        \right\}y\right]

    .. math::
        s(T) = 0 \text{ for } T \le a_7
        
    .. math::
        z = \frac{T}{T+a_6} \cdot \frac{a_7 + a_6}{a_7}

    .. math::
        y = \frac{T-a_7}{T+a_6} \text{ for } T > a_7 \text{ otherwise } 0
        
    Parameters
    ----------
    T : float
        Temperature [K]
    a1-a7 : float
        Coefficients
    J : float, optional
        Integral offset
        
    Returns
    -------
    S-S(0) : float
        Difference in entropy from 0 K , [J/mol/K]
    
    Notes
    -----
    Analytical integral as provided in [1]_ and verified with numerical
    integration. 
    
    Examples
    --------
    >>> TRCCp_integral_over_T(300, 4.0, 124000, 245, 50.539, -49.469, 
    ... 220440000, 560, 78)
    213.80156219151888
    
    References
    ----------
    .. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
       in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
    '''
    # Possible optimizations: pre-cache as much as possible.
    # If this were replaced by a cache, much of this would not need to be computed.
    if T <= a7:
        y = 0.
    else:
        y = (T - a7)/(T + a6)
    a6_inv = 1.0/a6
    x3 = a7 + a6
    z = T*x3/(a7*(T + a6))
    if T <= a7:
        s = 0.
    else:
        a72 = a7*a7
        a62_inv = a6_inv*a6_inv
        a7_a6 = a7*a6_inv # a7/a6
        a7_a6_2 = a7_a6*a7_a6
        a7_a6_4 = a7_a6_2*a7_a6_2
        x1 = (a4*a72 - a5)*a62_inv # part of third, sum
        first = (a3 + ((a4*a72 - a5)*a62_inv)*a7_a6_4)*a7_a6_2*log(z)
        second = (a3 + a4)*log((T + a6)/(x3))
        third = 0.0
        y_pow = 1.0
        a7_a6_pow = a7_a6_2*a7_a6_4
        na7_a6_inv = -1.0/a7_a6

        for i in range(1, 8):
            y_pow = y_pow*y
            a7_a6_pow *= na7_a6_inv
            third += (x1*a7_a6_pow - a4)*y_pow/i
        fourth = -(a3*a6_inv*x3 + a5*y_pow/(7.0*y*a7*x3))*y
        s = first + second + third + fourth
    x2 = a2/T
    return R*(J + a0*log(T) + a1/(a2*a2)*(1. + x2)*exp(-x2) + s)
    

### Heat capacities of liquids

def Rowlinson_Poling(T, Tc, omega, Cpgm):
    r'''Calculate liquid constant-pressure heat capacity with the [1]_ CSP method.
    This equation is not terrible accurate.

    The heat capacity of a liquid is given by:

    .. math::
        \frac{Cp^{L} - Cp^{g}}{R} = 1.586 + \frac{0.49}{1-T_r} +
        \omega\left[ 4.2775 + \frac{6.3(1-T_r)^{1/3}}{T_r} + \frac{0.4355}{1-T_r}\right]
        
    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor for fluid, [-]
    Cpgm : float
        Constant-pressure gas heat capacity, [J/mol/K]
        
    Returns
    -------
    Cplm : float
        Liquid constant-pressure heat capacity, [J/mol/K]
        
    Notes
    -----
    Poling compared 212 substances, and found error at 298K larger than 10%
    for 18 of them, mostly associating. Of the other 194 compounds, AARD is 2.5%.
    
    Examples
    --------
    >>> Rowlinson_Poling(350.0, 435.5, 0.203, 91.21)
    143.80196224081436
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    one_minus_Tr = 1. - Tr
    Cplm = Cpgm+ R*(1.586 + 0.49/one_minus_Tr + omega*(4.2775 + 6.3*one_minus_Tr**(1/3.)/Tr + 0.4355/one_minus_Tr))
    return Cplm

def Rowlinson_Bondi(T, Tc, omega, Cpgm):
    r'''Calculate liquid constant-pressure heat capacity with the CSP method
    shown in [1]_.
    
    The heat capacity of a liquid is given by:
        
    .. math::
        \frac{Cp^L - Cp^{ig}}{R} = 1.45 + 0.45(1-T_r)^{-1} + 0.25\omega
        [17.11 + 25.2(1-T_r)^{1/3}T_r^{-1} + 1.742(1-T_r)^{-1}]
        
    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor for fluid, [-]
    Cpgm : float
        Constant-pressure gas heat capacity, [J/mol/K]
        
    Returns
    -------
    Cplm : float
        Liquid constant-pressure heat capacity, [J/mol/K]
        
    Notes
    -----
    Less accurate than `Rowlinson_Poling`.
    
    Examples
    --------
    >>> Rowlinson_Bondi(T=373.28, Tc=535.55, omega=0.323, Cpgm=119.342)
    175.3976263003074
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
       Berlin; New York:: Springer, 2010.
    .. [3] J.S. Rowlinson, Liquids and Liquid Mixtures, 2nd Ed.,
       Butterworth, London (1969).
    '''
    Tr = T/Tc
    one_minus_Tr = 1. - Tr
    Cplm = Cpgm + R*(1.45 + 0.45/(one_minus_Tr) + 0.25*omega*(17.11 + 25.2*(one_minus_Tr)**(1/3.)/Tr + 1.742/one_minus_Tr))
    return Cplm

def Dadgostar_Shaw_terms(similarity_variable):
    """
    Return terms for the computation of Dadgostar-Shaw heat capacity equation.

    Parameters
    ----------
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]

    Returns
    -------
    first : float
        First term.
    second : float
        Second term.
    third : float
        Third term.

    See Also
    --------
    Dadgostar_Shaw

    """
    a = similarity_variable
    a2 = a*a
    a11 = -0.3416
    a12 = 2.2671
    a21 = 0.1064
    a22 = -0.3874
    a31 = -9.8231E-05
    a32 = 4.182E-04
    # Didn't seem to improve the comparison; sum of errors on some
    # points included went from 65.5  to 286.
    # Author probably used more precision in their calculation.
    #    constant = 3*R*(theta/T)**2*exp(theta/T)/(exp(theta/T)-1)**2
    constant = 24.5
    return (constant * (a11*a + a12*a2),
            a21*a + a22*a2,
            a31*a + a32*a2)

def Dadgostar_Shaw(T, similarity_variable, MW=None, terms=None):
    r'''Calculate liquid constant-pressure heat capacity with the similarity
    variable concept and method as shown in [1]_.
    
    .. math::
        C_{p} = 24.5(a_{11}\alpha + a_{12}\alpha^2)+ (a_{21}\alpha
        + a_{22}\alpha^2)T +(a_{31}\alpha + a_{32}\alpha^2)T^2
        
    Parameters
    ----------
    T : float
        Temperature of liquid [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]
    MW : float, optional
        Molecular weight of the pure compound or mixture average, [g/mol]
    terms : float, optional
        Terms in Dadgostar-Shaw equation as computed by :obj:`Dadgostar_Shaw_terms`
        
    Returns
    -------
    Cpl : float
        Liquid constant-pressure heat capacity, [J/mol/K if MW given; J/kg/K otherwise]
        
    Notes
    -----
    Many restrictions on its use.
    Original model is in terms of J/g/K. Note that the model is for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    a11 = -0.3416; a12 = 2.2671; a21 = 0.1064; a22 = -0.3874l;
    a31 = -9.8231E-05; a32 = 4.182E-04
    
    Examples
    --------
    >>> Dadgostar_Shaw(355.6, 0.139)
    1802.5291501191516
    
    References
    ----------
    .. [1] Dadgostar, Nafiseh, and John M. Shaw. "A Predictive Correlation for
       the Constant-Pressure Specific Heat Capacity of Pure and Ill-Defined
       Liquid Hydrocarbons." Fluid Phase Equilibria 313 (January 15, 2012):
       211-226. doi:10.1016/j.fluid.2011.09.015.
    '''
    first, second, third = terms or Dadgostar_Shaw_terms(similarity_variable)
    Cp = (first + second*T + third*T**2)
    return Cp*1000. if MW is None else Cp*MW

def Dadgostar_Shaw_integral(T, similarity_variable, MW=None, terms=None):
    r'''Calculate the integral of liquid constant-pressure heat capacity 
    with the similarity variable concept and method as shown in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature of gas [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]
    MW : float, optional
        Molecular weight of the pure compound or mixture average, [g/mol]
    terms : float, optional
        Terms in Dadgostar-Shaw equation as computed by :obj:`Dadgostar_Shaw_terms`
        
    Returns
    -------
    H : float
        Difference in enthalpy from 0 K, [J/mol if MW given; J/kg otherwise]
        
    Notes
    -----
    Original model is in terms of J/g/K. Note that the model is for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    Integral was computed with SymPy.
    
    See Also
    --------
    Dadgostar_Shaw
    Dadgostar_Shaw_integral_over_T
    
    Examples
    --------
    >>> Dadgostar_Shaw_integral(300.0, 0.1333)
    238908.15142664989
    
    References
    ----------
    .. [1] Dadgostar, Nafiseh, and John M. Shaw. "A Predictive Correlation for
       the Constant-Pressure Specific Heat Capacity of Pure and Ill-Defined
       Liquid Hydrocarbons." Fluid Phase Equilibria 313 (January 15, 2012):
       211-226. doi:10.1016/j.fluid.2011.09.015.
    
    '''
    T2 = T*T
    first, second, third = terms or Dadgostar_Shaw_terms(similarity_variable)
    H = T2*T/3.*third + T2*0.5*second + T*first
    return H*1000. if MW is None else H*MW

def Dadgostar_Shaw_integral_over_T(T, similarity_variable, MW=None, terms=None):
    r'''Calculate the integral of liquid constant-pressure heat capacity 
    with the similarity variable concept and method as shown in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature of gas [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]
    MW : float, optional
        Molecular weight of the pure compound or mixture average, [g/mol]
    terms : float, optional
        Terms in Dadgostar-Shaw equation as computed by :obj:`Dadgostar_Shaw_terms`
        
    Returns
    -------
    S : float
        Difference in entropy from 0 K, [J/mol/K if MW given; J/kg/K otherwise]
        
    Notes
    -----
    Original model is in terms of J/g/K. Note that the model is for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    Integral was computed with SymPy.
    
    See Also
    --------
    Dadgostar_Shaw
    Dadgostar_Shaw_integral
    
    Examples
    --------
    >>> Dadgostar_Shaw_integral_over_T(300.0, 0.1333)
    1201.1409113147918
    
    References
    ----------
    .. [1] Dadgostar, Nafiseh, and John M. Shaw. "A Predictive Correlation for
       the Constant-Pressure Specific Heat Capacity of Pure and Ill-Defined
       Liquid Hydrocarbons." Fluid Phase Equilibria 313 (January 15, 2012):
       211-226. doi:10.1016/j.fluid.2011.09.015.
    
    '''
    first, second, third = terms or Dadgostar_Shaw_terms(similarity_variable)
    S = T*T*0.5*third + T*second + first*log(T)
    return S*1000. if MW is None else S*MW

def Zabransky_quasi_polynomial(T, Tc, a1, a2, a3, a4, a5, a6):
    r'''Calculates liquid heat capacity using the model developed in [1]_.
    
    .. math::
        \frac{C}{R}=A_1\ln(1-T_r) + \frac{A_2}{1-T_r}
        + \sum_{j=0}^m A_{j+3} T_r^j
        
    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature of fluid, [K]
    a1-a6 : float
        Coefficients
        
    Returns
    -------
    Cp : float
        Liquid heat capacity, [J/mol/K]
        
    Notes
    -----
    Used only for isobaric heat capacities, not saturation heat capacities.
    Designed for reasonable extrapolation behavior caused by using the reduced
    critical temperature. Used by the authors of [1]_ when critical temperature
    was available for the fluid.
    Analytical integrals are available for this expression.
    
    Examples
    --------
    >>> Zabransky_quasi_polynomial(330, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    165.472878778683
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    Tr = T/Tc
    return R*(a1*log(1.0-Tr) + a2/(1.0-Tr) + a3 + Tr*(Tr*(Tr*a6 + a5) + a4))


def Zabransky_quasi_polynomial_integral(T, Tc, a1, a2, a3, a4, a5, a6):
    r'''Calculates the integral of liquid heat capacity using the  
    quasi-polynomial model developed in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature [K]
    a1-a6 : float
        Coefficients
        
    Returns
    -------
    H : float
        Difference in enthalpy from 0 K, [J/mol]
        
    Notes
    -----
    The analytical integral was derived with SymPy; it is a simple polynomial
    plus some logarithms.
    
    Examples
    --------
    >>> H2 = Zabransky_quasi_polynomial_integral(300, 591.79, -3.12743, 
    ... 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    >>> H1 = Zabransky_quasi_polynomial_integral(200, 591.79, -3.12743, 
    ... 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    >>> H2 - H1
    14662.031376528757
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    Tc2 = Tc*Tc
    Tc3 = Tc2*Tc
    term = T - Tc
    return R*(T*(T*(T*(T*a6/(4.*Tc3) + a5/(3.*Tc2)) + a4/(2.*Tc)) - a1 + a3) 
              + T*a1*log(1. - T/Tc) - 0.5*Tc*(a1 + a2)*log(term*term))

def Zabransky_quasi_polynomial_integral_over_T(T, Tc, a1, a2, a3, a4, a5, a6):
    r'''Calculates the integral of liquid heat capacity over T using the 
    quasi-polynomial model  developed in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature [K]
    a1-a6 : float
        Coefficients
        
    Returns
    -------
    S : float
        Difference in entropy from 0 K, [J/mol/K]
        
    Notes
    -----
    The analytical integral was derived with Sympy. It requires the 
    Polylog(2,x) function, which is unimplemented in SciPy. A very accurate 
    numerical approximation was implemented as :obj:`fluids.numerics.polylog2`.
    Relatively slow due to the use of that special function.
    
    Examples
    --------
    >>> S2 = Zabransky_quasi_polynomial_integral_over_T(300, 591.79, -3.12743, 
    ... 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    >>> S1 = Zabransky_quasi_polynomial_integral_over_T(200, 591.79, -3.12743, 
    ... 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    >>> S2 - S1
    59.16999297436473
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    term = T - Tc
    logT = log(T)
    Tc2 = Tc*Tc
    Tc3 = Tc2*Tc
    return R*(a3*logT -a1*polylog2(T/Tc) - a2*(-logT + 0.5*log(term*term))
              + T*(T*(T*a6/(3.*Tc3) + a5/(2.*Tc2)) + a4/Tc))

def Zabransky_cubic(T, a1, a2, a3, a4):
    r'''Calculates liquid heat capacity using the model developed in [1]_.
    
    .. math::
        \frac{C}{R}=\sum_{j=0}^3 A_{j+1} \left(\frac{T}{100 \text{K}}\right)^j
        
    Parameters
    ----------
    T : float
        Temperature [K]
    a1 : float
        Coefficient, [-]
    a2 : float
        Coefficient, [-]
    a3 : float
        Coefficient, [-]
    a4 : float
        Coefficient, [-]
        
    Returns
    -------
    Cp : float
        Liquid heat capacity, [J/mol/K]
        
    Notes
    -----
    Most often form used in [1]_.
    Analytical integrals are available for this expression.
    
    Examples
    --------
    >>> Zabransky_cubic(298.15, 20.9634, -10.1344, 2.8253, -0.256738)
    75.31465144297991
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    T = T/100.
    return R*(((a4*T + a3)*T + a2)*T + a1)

def Zabransky_cubic_integral(T, a1, a2, a3, a4):
    r'''Calculates the integral of liquid heat capacity using the model 
    developed in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature [K]
    a1 : float
        Coefficient, [-]
    a2 : float
        Coefficient, [-]
    a3 : float
        Coefficient, [-]
    a4 : float
        Coefficient, [-]
        
    Returns
    -------
    H : float
        Difference in enthalpy from 0 K, [J/mol]
        
    Notes
    -----
    The analytical integral was derived with Sympy; it is a simple polynomial.
    
    Examples
    --------
    >>> Zabransky_cubic_integral(298.15, 20.9634, -10.1344, 2.8253, -0.256738)
    31051.690370364562
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    T = T/100.
    return 100.0*R*T*(T*(T*(T*a4*0.25 + a3*(1.0/3.)) + a2*0.5) + a1)

def Zabransky_cubic_integral_over_T(T, a1, a2, a3, a4):
    r'''Calculates the integral of liquid heat capacity over T using the model 
    developed in [1]_.
    
    Parameters
    ----------
    T : float
        Temperature [K]
    a1 : float
        Coefficient, [-]
    a2 : float
        Coefficient, [-]
    a3 : float
        Coefficient, [-]
    a4 : float
        Coefficient, [-]
        
    Returns
    -------
    S : float
        Difference in entropy from 0 K, [J/mol/K]
        
    Notes
    -----
    The analytical integral was derived with Sympy; it is a simple polynomial,
    plus a logarithm
    
    Examples
    --------
    >>> Zabransky_cubic_integral_over_T(298.15, 20.9634, -10.1344, 2.8253, 
    ... -0.256738)
    24.732465342840914
    
    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    T = T*1e-2
    return R*(T*(T*(T*a4/3.0 + 0.5*a3) + a2) + a1*log(T))

### Solid

def Perry_151(T, a, b, c, d):
    r"""
    Return the solid molar heat capacity of a chemical using the Perry 151 method,
    as described in [1]_.
    
    Parameters
    ----------
    a,b,c,d : float
        Regressed coefficients.
    
    Returns
    -------
    Cps : float
        Solid constant-pressure heat capacity, [J/mol/K]
    
    Notes
    -----
    The solid heat capacity is given by:
    
    .. math:: C_n = 4.184 (a + bT + \frac{c}{T^2} + dT^2)
    
    Coefficients are listed in section 2, table 151 of [1]_. Note that the
    original model was in a Calorie basis, but has been translated to Joules.
    
    Examples
    --------
    Heat capacity of solid aluminum at 300 K:
        
    >>> Perry_151(300, 4.8, 0.00322, 0., 0.)
    24.124944
    
    References
    ----------
    .. [1] Green, Don, and Robert Perry.
       Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    
    """
    T2 = T**2
    return (a + b*T + c/T2 + d*T2) * 4.184

def Lastovka_solid(T, similarity_variable, MW=None):
    r'''Calculate solid constant-pressure heat capacity with the similarity
    variable concept and method as shown in [1]_.
    
    .. math::
        C_p = 3(A_1\alpha + A_2\alpha^2)R\left(\frac{\theta}{T}\right)^2
        \frac{\exp(\theta/T)}{[\exp(\theta/T)-1]^2}
        + (C_1\alpha + C_2\alpha^2)T + (D_1\alpha + D_2\alpha^2)T^2
        
    Parameters
    ----------
    T : float
        Temperature of solid [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]
    MW : float, optional
        Molecular weight of the pure compound or mixture average, [g/mol]
        
    Returns
    -------
    Cps : float
        Solid constant-pressure heat capacity, [J/mol/K if MW given; J/kg/K otherwise]
    
    Notes
    -----
    Many restrictions on its use. Trained on data with MW from 12.24 g/mol
    to 402.4 g/mol, C mass fractions from 61.3% to 95.2%,
    H mass fractions from 3.73% to 15.2%, N mass fractions from 0 to 15.4%,
    O mass fractions from 0 to 18.8%, and S mass fractions from 0 to 29.6%.
    Recommended for organic compounds with low mass fractions of hetero-atoms
    and especially when molar mass exceeds 200 g/mol. This model does not show
    and effects of phase transition but should not be used passed the triple
    point.
    Original model is in terms of J/g/K. Note that the model s for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    
    A1 = 0.013183
    
    A2 = 0.249381
    
    :math:`\theta` = 151.8675
    
    C1 = 0.026526
    
    C2 = -0.024942
    
    D1 = 0.000025
    
    D2 = -0.000123
    
    Examples
    --------
    >>> Lastovka_solid(300, 0.2139)
    1682.0637469909211
    
    References
    ----------
    .. [1] Laštovka, Václav, Michal Fulem, Mildred Becerra, and John M. Shaw.
       "A Similarity Variable for Estimating the Heat Capacity of Solid Organic
       Compounds: Part II. Application: Heat Capacity Calculation for
       Ill-Defined Organic Solids." Fluid Phase Equilibria 268, no. 1-2
       (June 25, 2008): 134-41. doi:10.1016/j.fluid.2008.03.018.
    '''
    A1 = 0.013183
    A2 = 0.249381
    theta = 151.8675
    C1 = 0.026526
    C2 = -0.024942
    D1 = 0.000025
    D2 = -0.000123
    theta_div_T = theta/T
    exp_term = exp(theta_div_T)
    a = similarity_variable
    Cp = a*(3.0*(A1 + A2*a)*R*(theta_div_T)**2*exp_term/(exp_term-1)**2
          + (C1 + C2*a)*T
          + (D1 + D2*a)*T**2)
    return Cp*1000. if MW is None else Cp*MW

def Lastovka_solid_integral(T, similarity_variable, MW=None):
    r'''Integrates solid constant-pressure heat capacity with the similarity
    variable concept and method as shown in [1]_.
    
    uses an explicit form as derived with Sympy.
    
    Parameters
    ----------
    T : float
        Temperature of solid [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]
    MW : float, optional
        Molecular weight of the pure compound or mixture average, [g/mol]
        
    Returns
    -------
    H : float
        Difference in enthalpy from 0 K, [J/mol if MW given; J/kg otherwise]
        
    Notes
    -----
    Original model is in terms of J/g/K. Note that the model is for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    
    See Also
    --------
    Lastovka_solid
    
    Examples
    --------
    >>> Lastovka_solid_integral(300, 0.2139)
    283246.1519409122
    
    References
    ----------
    .. [1] Laštovka, Václav, Michal Fulem, Mildred Becerra, and John M. Shaw.
       "A Similarity Variable for Estimating the Heat Capacity of Solid Organic
       Compounds: Part II. Application: Heat Capacity Calculation for
       Ill-Defined Organic Solids." Fluid Phase Equilibria 268, no. 1-2
       (June 25, 2008): 134-41. doi:10.1016/j.fluid.2008.03.018.
    '''
    A1 = 0.013183
    A2 = 0.249381
    theta = 151.8675
    C1 = 0.026526
    C2 = -0.024942
    D1 = 0.000025
    D2 = -0.000123
    a = similarity_variable
    T2 = T*T
    H = a*(T*T2*(D1 + D2*a)/3.
           + 0.5*T2*(C1 + C2*a)
           + 3.0*R*theta*(A1 + A2*a)/(exp(theta/T) - 1.))
    return H*1000. if MW is None else H*MW

def Lastovka_solid_integral_over_T(T, similarity_variable, MW=None):
    r'''Integrates over T solid constant-pressure heat capacity with the 
    similarity variable concept and method as shown in [1]_.
    
    uses an explicit form as derived with Sympy.
    
    Parameters
    ----------
    T : float
        Temperature of solid [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]
    MW : float, optional
        Molecular weight of the pure compound or mixture average, [g/mol]
        
    Returns
    -------
    S : float
        Difference in entropy from 0 K, [J/mol/K if MW given; J/kg/K otherwise]
    
    Notes
    -----
    Original model is in terms of J/g/K. Note that the model is for predicting
    mass heat capacity, not molar heat capacity like most other methods!
    
    See Also
    --------
    Lastovka_solid
    
    Examples
    --------
    >>> Lastovka_solid_integral_over_T(300, 0.2139)
    1947.5537561495564
    
    References
    ----------
    .. [1] Laštovka, Václav, Michal Fulem, Mildred Becerra, and John M. Shaw.
       "A Similarity Variable for Estimating the Heat Capacity of Solid Organic
       Compounds: Part II. Application: Heat Capacity Calculation for
       Ill-Defined Organic Solids." Fluid Phase Equilibria 268, no. 1-2
       (June 25, 2008): 134-41. doi:10.1016/j.fluid.2008.03.018.
    '''
    A1 = 0.013183
    A2 = 0.249381
    theta = 151.8675
    C1 = 0.026526
    C2 = -0.024942
    D1 = 0.000025
    D2 = -0.000123
    
    a = similarity_variable
    exp_theta_T = exp(theta/T)
    A_term = (A1 + A2*a)
    S = a*(-3.0*R*A_term*log(exp_theta_T - 1.) 
           + 0.5*T**2*(D1 + D2*a)
           + T*(C1 + C2*a)
           + 3.0*R*theta*A_term*(1/(T*exp_theta_T - T) + 1/T))
    return S*1000. if MW is None else S*MW
