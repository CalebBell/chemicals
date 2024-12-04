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

This module contains miscellaneous functions which may be useful. This includes
definitions of some chemical properties, and conversions between others.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

"""

__all__ = ['isobaric_expansion', 'isothermal_compressibility',
'Cp_minus_Cv', 'speed_of_sound', 'Joule_Thomson',
'phase_identification_parameter', 'phase_identification_parameter_phase',
'isentropic_exponent', 'isentropic_exponent_TV', 'isentropic_exponent_PT', 'isentropic_exponent_PV',
'Vm_to_rho', 'rho_to_Vm',
'Z',  'zs_to_ws', 'ws_to_zs', 'zs_to_Vfs',
'Vfs_to_zs',
'ms_to_ns', 'ns_to_ms', 'ns_to_Qls', 'Qls_to_ns',
'Qls_to_ms', 'ms_to_Qls',
 'none_and_length_check', 'normalize', 'remove_zeros',
 'mixing_simple',
'mixing_logarithmic', 'mixing_power', 'to_num', 'Parachor', 'property_molar_to_mass', 'property_mass_to_molar',
'SG_to_API', 'API_to_SG', 'API_to_rho', 'rho_to_API', 'SG',   'Watson_K',
'dxs_to_dns', 'dns_to_dn_partials', 'dxs_to_dn_partials', 'd2ns_to_dn2_partials',
'd2xs_to_dxdn_partials', 'dxs_to_dxsn1', 'd2xs_to_d2xsn1',
 'vapor_mass_quality', 'mix_component_flows',
'mix_multiple_component_flows', 'mix_component_partial_flows',
'solve_flow_composition_mix', 'radius_of_gyration',
'v_to_v_molar', 'v_molar_to_v', 'molar_velocity_to_velocity', 'velocity_to_molar_velocity']

import os
import sys
from math import (  # Not supported in Python 2.6: expm1, erf, erfc,gamma lgamma
    exp,
    pi,
    sqrt,
)

# __all__.extend(['R', 'k', 'N_A', 'calorie', 'epsilon_0']) # 'expm1', 'erf', 'erfc',  'lgamma', 'gamma',
# Obtained from SciPy 0.19 (2014 CODATA)
# Included here so calculations are consistent across SciPy versions
from fluids.constants import N_A, R
from fluids.numerics import numpy as np
from fluids.numerics import trunc_log, trunc_exp, cbrt

__all__.extend(['PY37'])
version_components = sys.version.split('.')
PY_MAJOR, PY_MINOR = int(version_components[0]), int(version_components[1])
PY37 = (PY_MAJOR, PY_MINOR) >= (3, 7)

try:
    source_path = os.path.dirname(__file__) # micropython
except:
    source_path = ''

if os.name == 'nt':
    def os_path_join(*args):
        return '\\'.join(args)
else:
    def os_path_join(*args):
        return '/'.join(args)

can_load_data = True
try:
    implementation = sys.implementation.name
    if implementation  in ('micropython', 'ironpython'):
        can_load_data = False
except:
    pass

numba_blacklisted = ['mark_numba_incompatible', 'mark_numba_uncacheable']
numba_cache_blacklisted = []

def mark_numba_incompatible(f):
    numba_blacklisted.append(f.__name__)
    return f

def mark_numba_uncacheable(f):
    numba_cache_blacklisted.append(f.__name__)
    return f

@mark_numba_incompatible
def to_num(values):
    r'''Legacy function to turn a list of strings into either floats
    (if numeric), stripped strings (if not) or None if the string is empty.
    Accepts any numeric formatting the float function does.

    Parameters
    ----------
    values : list
        list of strings

    Returns
    -------
    values : list
        list of floats, strings, and None values [-]

    Examples
    --------
    >>> to_num(['1', '1.1', '1E5', '0xB4', ''])
    [1.0, 1.1, 100000.0, '0xB4', None]
    '''
    float_ = float
    for i in range(len(values)):
        try:
            values[i] = float_(values[i])
        except:
            if values[i] == '':
                values[i] = None
            else:
                values[i] = values[i].strip()
    return values


try:
    ndarray = np.ndarray
except:
    from fluids.numerics import FakePackage as ndarray




empty_dict = {}
immutable_types = {type(None), bool, int, float, complex, str, bytes,type,
                       # dictionary views are read-only, kinda unexpected
                      type(empty_dict.keys()), type(empty_dict.items()), type(empty_dict.values()),
                      }
try:
    numpy_immutable_types = [np.bool_, np.byte, np.ubyte, np.short, np.ushort, np.intc, np.uintc, np.int_, np.uint, np.longlong, np.ulonglong, np.half, np.float16, np.single, np.double, np.longdouble, np.csingle, np.cdouble, np.clongdouble]
    immutable_types.update(numpy_immutable_types)
except:
    pass

immutable_class_types = {'Fraction', 'Decimal'}

object_data_type_cache = {}
def create_object_data_function(instance):
    t = type(instance)
    try:
        return object_data_type_cache[t]
    except:
        pass
    has_d = hasattr(instance, '__dict__')
    slots = []
    for base in instance.__class__.__mro__[:-1]:
        try:
            slots.extend(base.__slots__)
        except:
            continue
    # remove duplicates
    slots_set = set(slots)
    try:
        slots_set.remove('__dict__')
    except:
        pass
    slots = tuple(slots_set)

    if has_d and slots:
        def specialized_object_data(obj, slots=slots):
            d = obj.__dict__.copy()
            for s in slots:
                try:
                    d[s] = getattr(obj, s)
                except:
                    pass
            return d
    elif has_d:
        def specialized_object_data(obj):
            return obj.__dict__.copy()
    else:
        def specialized_object_data(obj, slots=slots):
            d = {}
            for s in slots:
                try:
                    d[s] = getattr(obj, s)
                except:
                    pass
            return d

    object_data_type_cache[t] = specialized_object_data
    return specialized_object_data

@mark_numba_incompatible
def object_data(obj):
    t = type(obj)
    try:
        return object_data_type_cache[t](obj)
    except:
        return create_object_data_function(obj)(obj)
    # Always returns a new dictionary
    # try:
    #     d = obj.__dict__.copy()
    # except:
    #     d = {}

    # for base in obj.__class__.__mro__[:-1]:
    #     try:
    #         slots = base.__slots__
    #     except:
    #         continue
    #     for s in slots:
    #         if s == '__dict__':
    #             continue
    #         try:
    #             d[s] = getattr(obj, s)
    #         except:
    #             pass
    # return d

@mark_numba_incompatible
def recursive_copy(obj):
    obj_type = type(obj)
    if obj_type in immutable_types:
        return obj
    elif obj_type is tuple:
        return tuple(recursive_copy(v) for v in obj)
    elif obj_type is list:
        return [recursive_copy(v) for v in obj]
    elif obj_type is dict:
        d = {}
        for k, v in obj.items():
            d[recursive_copy(k)] = recursive_copy(v)
        return d
    elif obj_type is bytearray:
        return obj.copy()
    elif obj_type is np.ndarray:
        # Mutable objects
        # .copy() won't work for obj arrays of mutable values, so have to handle it separately
        if obj.dtype == object:
            copy = np.array([recursive_copy(o) for o in obj.ravel().tolist()], dtype=object)
            copy.reshape(obj.shape)
            return copy
        return obj.copy()
    elif obj_type is range:
        return range(obj.start, obj.stop, obj.step)
    elif obj_type is set:
        return {recursive_copy(v) for v in obj}
    elif obj_type is frozenset:
        return frozenset([recursive_copy(v) for v in obj])

    # We need to handle the copy on some object types
    # that we don't want to import, so use their class names as a good proxy

    obj_class_name = obj_type.__name__

    if obj_class_name == 'array':
        return obj.__copy__() # mutable
    elif obj_class_name in immutable_class_types:
        return obj

    if hasattr(obj, '__copy__'):
        # Allow objects to provide their own hooks and wash our hands of them
        return obj.__copy__()

    # Is it a named tuple
    if hasattr(obj, '_fields') and hasattr(obj, '_asdict') and obj_type.__mro__[1] is tuple:
        return obj_type(*(recursive_copy(v) for v in obj))
    raise ValueError(f"No copy function implemented for type {obj_type}")

@mark_numba_incompatible
def hash_any_primitive(v):
    '''Method to hash a primitive - with basic support for lists and
    dictionaries.

    Parameters
    ----------
    v : object
        Value to hash [-]

    Returns
    -------
    h : int
        Hashed value [-]

    Notes
    -----
    Handles up to 3d lists. Assumes all lists are the same dimension.
    Handles dictionaries by iterating over each value and key.

    Will fail explosively on circular references.

    The values returned by this function should not be counted on to be the
    same in the future.

    Values change on startup due to hash randomization on python 3.3+.

    Examples
    --------
    hash_any_primitive([1,2,3,4,5])

    hash_any_primitive({'a': [1,2,3], 'b': []})
    '''
    t = type(v)
    if t is list:
        if len(v) and isinstance(v[0], list):
            if len(v[0]) and isinstance(v[0][0], list):
                # 3d
                v = tuple(tuple(hash_any_primitive(j) for j in i) for i in v)
            else:
                # 2d
                v = tuple(hash_any_primitive(i) for i in v)
        else:
            # likely a 1d list
            v = tuple(i if type(i) in immutable_types else hash_any_primitive(i) for i in v)
    elif t is dict:
        temp_hash = []
        # Do not want to order this at all
        for key, value in v.items():
            # Do not bother hashing the value if it's immutable as it will be hashed with the tuple
            value_hash = value if type(value) in immutable_types else hash_any_primitive(value)
            key_value_hash = hash((key, value_hash))
            temp_hash.append(key_value_hash)
        v = hash(frozenset(temp_hash))
    elif t is set:
        # Should only contain hashable items
        v = frozenset(v)
    elif t is tuple:
        v = tuple(i if type(i) in immutable_types else hash_any_primitive(i) for i in v)
    elif t is ndarray:
        v = hash(v.data.tobytes())
    return hash(v)

def Parachor(MW, rhol, rhog, sigma):
    r'''Calculate Parachor for a pure species, using its density in the
    liquid and gas phases, surface tension, and molecular weight.

    .. math::
        P = \frac{\sigma^{0.25} MW}{\rho_L - \rho_V}

    Parameters
    ----------
    MW : float
        Molecular weight, [g/mol]
    rhol : float
        Liquid density [kg/m^3]
    rhog : float
        Gas density [kg/m^3]
    sigma : float
        Surface tension, [N/m]

    Returns
    -------
    P : float
        Parachor, [N^0.25*m^2.75/mol]

    Notes
    -----
    To convert the output of this function to units of [mN^0.25*m^2.75/kmol],
    multiply by 5623.4132519.

    Values in group contribution tables for Parachor are often listed as
    dimensionless, in which they are multiplied by 5623413 and the appropriate
    units to make them dimensionless.

    Examples
    --------
    Calculating Parachor from a known surface tension for methyl isobutyl
    ketone at 293.15 K

    >>> Parachor(100.15888, 800.8088185536124, 4.97865317223119, 0.02672166960656005)
    5.088443542210164e-05

    Converting to the `dimensionless` form:

    >>> 5623413*5.088443542210164e-05
    286.14419565030687

    Compared to 274.9 according to a group contribution method described in
    [3]_.

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       8E. McGraw-Hill Professional, 2007.
    .. [3] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    rhol, rhog = rhol*1000., rhog*1000. # Convert kg/m^3 to g/m^3
    return sigma**0.25*MW/(rhol-rhog) # (N/m)**0.25*g/mol/(g/m^3)


def property_molar_to_mass(A_molar, MW):
    r'''Convert a quantity in molar units [thing/mol] to mass units [thing/kg].
    The standard gram-mole is used here, as it is everwhere in this library.

    .. math::
        A_{\text{mass}} = \frac{1000 A_{\text{molar}}}{\text{MW}}

    Parameters
    ----------
    A_molar : float
        Quantity in molar units [thing/mol]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    A_mass : float
        Quantity in molar units [thing/kg]

    Notes
    -----
    For legacy reasons, if the value `A_molar` is None, None is also returned
    and no exception is returned.

    Examples
    --------
    >>> property_molar_to_mass(500, 18.015)
    27754.648903691366
    '''
    if A_molar is None:
        return None
    return A_molar*1000.0/MW


def property_mass_to_molar(A_mass, MW):
    r'''Convert a quantity in mass units [thing/kg] to molar units [thing/mol].
    The standard gram-mole is used here, as it is everwhere in this library.

    .. math::
        A_{\text{molar}} = \frac{A_{\text{mass}} \text{MW}}{1000}

    Parameters
    ----------
    A_mass : float
        Quantity in molar units [thing/kg]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    A_molar : float
        Quantity in molar units [thing/mol]

    Notes
    -----
    For legacy reasons, if the value `A_mass` is None, None is also returned
    and no exception is returned.

    Examples
    --------
    >>> property_mass_to_molar(20.0, 18.015)
    0.3603
    '''
    if A_mass is None:
        return None
    return 1e-3*A_mass*MW

root_1000 = 1000**0.5
root_1000_inv = 1.0/root_1000

def v_to_v_molar(v, MW):
    r'''Convert a velocity from units of m/s to a "molar" form of velocity,
    compatible with thermodynamic calculations on a molar basis.

    .. math::
        v\left(\frac{\text{m}\sqrt{\text{kg}} }{s \sqrt{\text{mol}}} \right)
        = v \text{(m/s)}
        \sqrt{\text{MW (g/mol)}}\cdot
        \left(\frac{1000 \text{g}}{1 \text{kg}}\right)^{-0.5}

    Parameters
    ----------
    v : float
        Velocity, [m/s]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    v_molar : float
        Molar velocity, [m*kg^0.5/s/mol^0.5]

    Examples
    --------
    >>> v_to_v_molar(500, 18.015)
    67.10998435404377
    '''
    return v*MW**0.5*root_1000_inv

def v_molar_to_v(v_molar, MW):
    r'''Convert a velocity from units of the molar velocity form to standard
    m/s units.

    .. math::
        v \text{(m/s)} = v\left(\frac{\text{m}\sqrt{\text{kg}} }
        {s \sqrt{\text{mol}}} \right)
        {\text{MW (g/mol)}}^{-0.5}\cdot
        \left(\frac{1000 \text{g}}{1 \text{kg}}\right)^{0.5}

    Parameters
    ----------
    v_molar : float
        Molar velocity, [m*kg^0.5/s/mol^0.5]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    v : float
        Velocity, [m/s]

    Examples
    --------
    >>> v_molar_to_v(67.10998435404377, 18.015)
    499.99999999999994
    '''
    return v_molar*root_1000*MW**-0.5

def vapor_mass_quality(VF, MWl, MWg):
    r'''Calculates the vapor quality on a mass basis of a two-phase mixture;
    this is the most common definition, where 1 means a pure vapor and 0 means
    a pure liquid. The vapor quality on a mass basis is related to the mole
    basis vapor fraction according to the following relationship:

    .. math::
        x = \frac{\frac{V}{F}\cdot  \text{MW}_g}
        {(1-\frac{V}{F})\text{MW}_l + \frac{V}{F}\text{MW}_g}

    Parameters
    ----------
    VF : float
        Mole-basis vapor fraction (0 = pure vapor, 1 = pure liquid), [-]
    MWl : float
        Average molecular weight of the liquid phase, [g/mol]
    MWg : float
        Average molecular weight of the vapor phase, [g/mol]

    Returns
    -------
    quality : float
        Vapor mass fraction of the two-phase system, [-]

    Notes
    -----
    Other definitions of vapor fraction use an enthalpy basis instead of a mass
    basis; still other less common ones take 1 to be the value of the
    liquid, and 0 as pure vapor.

    Examples
    --------
    >>> vapor_mass_quality(0.5, 60, 30)
    0.3333333333333333

    References
    ----------
    .. [1] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       8E. McGraw-Hill Professional, 2007.
    '''
    ng = VF
    nl = (1. - VF)
    return ng*MWg/(nl*MWl + ng*MWg)

def rho_to_API(rho, rho_ref=999.0170824078306):
    r'''Calculates API of a liquid given its mass density, as shown in
    [1]_.

    .. math::
        \text{API gravity} = \frac{141.5\rho_{ref}}{\rho} - 131.5

    Parameters
    ----------
    rho : float
        Mass density the fluid at 60 degrees Farenheight [kg/m^3]
    rho_ref : float, optional
        Density of the reference substance, [kg/m^3]

    Returns
    -------
    API : float
        API of the fluid [-]

    Notes
    -----
    Defined only at 60 degrees Fahrenheit.

    Examples
    --------
    >>> rho_to_API(820)
    40.8913623
    >>> SG_to_API(SG(820))
    40.8913623

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    '''
    return 141.5*rho_ref/rho - 131.5

def API_to_rho(API, rho_ref=999.0170824078306):
    r'''Calculates mass density of a liquid given its API, as shown in
    [1]_.

    .. math::
        \rho~60^\circ\text{F} = \frac{141.5\rho_{ref}}{\text{API} + 131.5}

    Parameters
    ----------
    API : float
        API of the fluid [-]
    rho_ref : float, optional
        Density of the reference substance, [kg/m^3]

    Returns
    -------
    rho : float
        Mass density the fluid at 60 degrees Farenheight [kg/m^3]

    Notes
    -----
    Defined only at 60 degrees Fahrenheit.

    Examples
    --------
    >>> API_to_rho(rho_to_API(820))
    820.0

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    '''
    return (141.5*rho_ref)/(API + 131.5)



def SG_to_API(SG):
    r'''Calculates API of a liquid given its specific gravity, as shown in
    [1]_.

    .. math::
        \text{API gravity} = \frac{141.5}{\text{SG}} - 131.5

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]

    Returns
    -------
    API : float
        API of the fluid [-]

    Notes
    -----
    Defined only at 60 degrees Fahrenheit.

    Examples
    --------
    >>> SG_to_API(0.7365)
    60.62491513917175

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    '''
    return 141.5/SG - 131.5


def API_to_SG(API):
    r'''Calculates specific gravity of a liquid given its API, as shown in
    [1]_.

    .. math::
        \text{SG at}~60^\circ\text{F} = \frac{141.5}{\text{API gravity} +131.5}

    Parameters
    ----------
    API : float
        API of the fluid [-]

    Returns
    -------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]

    Notes
    -----
    Defined only at 60 degrees Fahrenheit.

    Examples
    --------
    >>> API_to_SG(60.62)
    0.7365188423901728

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    '''
    return 141.5/(API + 131.5)


def SG(rho, rho_ref=999.0170824078306):
    r'''Calculates the specific gravity of a substance with respect to another
    substance; by default, this is water at 15.555 °C (60 °F). For gases,
    normally the reference density is 1.2 kg/m^3, that of dry air. However, in
    general specific gravity should always be specified with respect to the
    temperature and pressure of its reference fluid. This can vary widely.

    .. math::
        SG = \frac{\rho}{\rho_{ref}}

    Parameters
    ----------
    rho : float
        Density of the substance, [kg/m^3]
    rho_ref : float, optional
        Density of the reference substance, [kg/m^3]

    Returns
    -------
    SG : float
        Specific gravity of the substance with respect to the reference
        density, [-]

    Notes
    -----
    Another common reference point is water at 4°C (rho_ref=999.9748691393087).
    Specific gravity is often used by consumers instead of density.
    The reference for solids is normally the same as for liquids - water.

    Examples
    --------
    >>> SG(860)
    0.8608461408159591
    '''
    return rho/rho_ref


def Watson_K(Tb, SG):
    r'''Calculates the Watson or UOP K Characterization factor
    of a liquid of a liquid given its specific gravity, and its
    average boiling point as shown in [1]_.

    .. math::
        K_W = \frac{T_b^{1/3}}{\text{SG at}~60^\circ\text{F}}

    Parameters
    ----------
    Tb : float
        Average normal boiling point, [K]
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]

    Returns
    -------
    K_W : float
        Watson characterization factor

    Notes
    -----
    There are different ways to compute the average boiling point,
    so two different definitions are often used - K_UOP using volume
    average boiling point (VABP) using distillation points of 10%, 30%,
    50%, 70%, and 90%; and K_Watson using mean average boiling point (MeABP).

    Examples
    --------
    >>> Watson_K(400, .8)
    11.20351186639291

    Sample problem in Comments on Procedure 2B5.1 of [1]_;
    a fluids has a MEAB of 580 F and a SG of 34.5.

    >>> from fluids.core import F2K
    >>> Watson_K(F2K(580), API_to_SG(34.5))
    11.884570347084471

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    '''
    return (Tb*1.8)**(1.0/3.0)/SG


def isobaric_expansion(V, dV_dT):
    r'''Calculate the isobaric coefficient of a thermal expansion, given its
    molar volume at a certain `T` and `P`, and its derivative of molar volume
    with respect to `T`.

    .. math::
        \beta = \frac{1}{V}\left(\frac{\partial V}{\partial T} \right)_P

    Parameters
    ----------
    V : float
        Molar volume at `T` and `P`, [m^3/mol]
    dV_dT : float
        Derivative of molar volume with respect to `T`, [m^3/mol/K]

    Returns
    -------
    beta : float
        Isobaric coefficient of a thermal expansion, [1/K]

    Notes
    -----
    For an ideal gas, this expression simplified to:

    .. math::
        \beta = \frac{1}{T}

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):

    >>> isobaric_expansion(0.000130229900873546, 1.58875261849113e-7)
    0.0012199599384121608

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return dV_dT/V


def isothermal_compressibility(V, dV_dP):
    r'''Calculate the isothermal coefficient of compressibility, given its
    molar volume at a certain `T` and `P`, and its derivative of molar volume
    with respect to `P`.

    .. math::
        \kappa = -\frac{1}{V}\left(\frac{\partial V}{\partial P} \right)_T

    Parameters
    ----------
    V : float
        Molar volume at `T` and `P`, [m^3/mol]
    dV_dP : float
        Derivative of molar volume with respect to `P`, [m^3/mol/Pa]

    Returns
    -------
    kappa : float
        Isothermal coefficient of compressibility, [1/Pa]

    Notes
    -----
    For an ideal gas, this expression simplified to:

    .. math::
        \kappa = \frac{1}{P}


    The isothermal bulk modulus is the inverse of this quantity:

    .. math::
        K = -V\left(\frac{\partial P}{\partial V} \right)_T

    The ideal gas isothermal bulk modulus is simply the gas's pressure.

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):

    >>> isothermal_compressibility(0.000130229900873546, -2.72902118209903e-13)
    2.095541165119158e-09

    Calculate the bulk modulus of propane from the PR EOS at 294 K as a gas:

    >>> 1/isothermal_compressibility(0.0024576770482135617, -3.5943321700795866e-09)
    683764.5859979445

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return -dV_dP/V


def phase_identification_parameter(V, dP_dT, dP_dV, d2P_dV2, d2P_dVdT):
    r'''Calculate the Phase Identification Parameter developed in [1]_ for
    the accurate and efficient determination of whether a fluid is a liquid or
    a gas based on the results of an equation of state. For supercritical
    conditions, this provides a good method for choosing which property
    correlations to use.

    .. math::
        \Pi = V \left[\frac{\frac{\partial^2 P}{\partial V \partial T}}
        {\frac{\partial P }{\partial T}}- \frac{\frac{\partial^2 P}{\partial
        V^2}}{\frac{\partial P}{\partial V}} \right]

    Parameters
    ----------
    V : float
        Molar volume at `T` and `P`, [m^3/mol]
    dP_dT : float
        Derivative of `P` with respect to `T`, [Pa/K]
    dP_dV : float
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]
    d2P_dV2 : float
        Second derivative of `P` with respect to `V`, [Pa*mol^2/m^6]
    d2P_dVdT : float
        Second derivative of `P` with respect to both `V` and `T`, [Pa*mol/m^3/K]

    Returns
    -------
    PIP : float
        Phase Identification Parameter, [-]

    Notes
    -----
    Heuristics were used by process simulators before the invent of this
    parameter.

    The criteria for liquid is Pi > 1; for vapor, Pi <= 1.

    There is also a solid phase mechanism available. For solids, the Solid
    Phase Identification Parameter is greater than 1, like liquids; however,
    unlike liquids, d2P_dVdT is always >0; it is < 0 for liquids and gases.

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):

    >>> phase_identification_parameter(0.000130229900874, 582169.397484,
    ... -3.66431747236e+12, 4.48067893805e+17, -20518995218.2)
    11.33428990564796

    References
    ----------
    .. [1] Venkatarathnam, G., and L. R. Oellrich. "Identification of the Phase
       of a Fluid Using Partial Derivatives of Pressure, Volume, and
       Temperature without Reference to Saturation Properties: Applications in
       Phase Equilibria Calculations." Fluid Phase Equilibria 301, no. 2
       (February 25, 2011): 225-33. doi:10.1016/j.fluid.2010.12.001.
    .. [2] Jayanti, Pranava Chaitanya, and G. Venkatarathnam. "Identification
       of the Phase of a Substance from the Derivatives of Pressure, Volume and
       Temperature, without Prior Knowledge of Saturation Properties: Extension
       to Solid Phase." Fluid Phase Equilibria 425 (October 15, 2016): 269-277.
       doi:10.1016/j.fluid.2016.06.001.
    '''
    return V*(d2P_dVdT/dP_dT - d2P_dV2/dP_dV)


def phase_identification_parameter_phase(d2P_dVdT, V=None, dP_dT=None, dP_dV=None, d2P_dV2=None):
    r'''Uses the Phase Identification Parameter concept developed in [1]_ and
    [2]_ to determine if a chemical is a solid, liquid, or vapor given the
    appropriate thermodynamic conditions.

    The criteria for liquid is PIP > 1; for vapor, PIP <= 1.

    For solids, PIP(solid) is defined to be d2P_dVdT. If it is larger than 0,
    the species is a solid. It is less than 0 for all liquids and gases.

    Parameters
    ----------
    d2P_dVdT : float
        Second derivative of `P` with respect to both `V` and `T`, [Pa*mol/m^3/K]
    V : float, optional
        Molar volume at `T` and `P`, [m^3/mol]
    dP_dT : float, optional
        Derivative of `P` with respect to `T`, [Pa/K]
    dP_dV : float, optional
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]
    d2P_dV2 : float, optional
        Second derivative of `P` with respect to `V`, [Pa*mol^2/m^6]

    Returns
    -------
    phase : str
        Either 's', 'l' or 'g'

    Notes
    -----
    The criteria for being a solid phase is checked first, which only
    requires d2P_dVdT. All other inputs are optional for this reason.
    However, an exception will be raised if the other inputs become
    needed to determine if a species is a liquid or a gas.

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):

    >>> phase_identification_parameter_phase(-20518995218.2, 0.000130229900874,
    ... 582169.397484, -3.66431747236e+12, 4.48067893805e+17)
    'l'

    References
    ----------
    .. [1] Venkatarathnam, G., and L. R. Oellrich. "Identification of the Phase
       of a Fluid Using Partial Derivatives of Pressure, Volume, and
       Temperature without Reference to Saturation Properties: Applications in
       Phase Equilibria Calculations." Fluid Phase Equilibria 301, no. 2
       (February 25, 2011): 225-33. doi:10.1016/j.fluid.2010.12.001.
    .. [2] Jayanti, Pranava Chaitanya, and G. Venkatarathnam. "Identification
       of the Phase of a Substance from the Derivatives of Pressure, Volume and
       Temperature, without Prior Knowledge of Saturation Properties: Extension
       to Solid Phase." Fluid Phase Equilibria 425 (October 15, 2016): 269-277.
       doi:10.1016/j.fluid.2016.06.001.
    '''
    if d2P_dVdT > 0:
        return 's'
    else:
        PIP = phase_identification_parameter(V=V, dP_dT=dP_dT, dP_dV=dP_dV,
                                             d2P_dV2=d2P_dV2, d2P_dVdT=d2P_dVdT)
        return 'l' if PIP > 1 else 'g'


def Cp_minus_Cv(T, dP_dT, dP_dV):
    r'''Calculate the difference between a real gas's constant-pressure heat
    capacity and constant-volume heat capacity, as given in [1]_, [2]_, and
    [3]_. The required derivatives should be calculated with an equation of
    state.

    .. math::
        C_p - C_v = -T\left(\frac{\partial P}{\partial T}\right)_V^2/
        \left(\frac{\partial P}{\partial V}\right)_T

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    dP_dT : float
        Derivative of `P` with respect to `T`, [Pa/K]
    dP_dV : float
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]

    Returns
    -------
    Cp_minus_Cv : float
        Cp - Cv for a real gas, [J/mol/K]

    Notes
    -----
    Equivalent expressions are:

    .. math::
        C_p - C_v= -T\left(\frac{\partial V}{\partial T}\right)_P^2/\left(
        \frac{\partial V}{\partial P}\right)_T

    .. math::
        C_p - C_v = T\left(\frac{\partial P}{\partial T}\right)
        \left(\frac{\partial V}{\partial T}\right)

    Note that these are not second derivatives, only first derivatives, some
    of which are squared.

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):

    >>> Cp_minus_Cv(299, 582232.475794113, -3665180614672.253)
    27.654681381642394

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering.
       Butterworth-Heinemann, 1985.
    .. [3] Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey.
       Chemical Thermodynamics for Process Simulation. 1st edition. Weinheim:
       Wiley-VCH, 2012.
    '''
    return -T*dP_dT*dP_dT/dP_dV


def speed_of_sound(V, dP_dV, Cp, Cv, MW=None):
    r'''Calculate a real fluid's speed of sound. The required derivatives should
    be calculated with an equation of state, and `Cp` and `Cv` are both the
    real fluid versions. Expression is given in [1]_ and [2]_; a unit conversion
    is further performed to obtain a result in m/s. If MW is not provided the
    result is returned in units of m*kg^0.5/s/mol^0.5.

    .. math::
        w = \left[-V^2 \left(\frac{\partial P}{\partial V}\right)_T \frac{C_p}
        {C_v}\right]^{1/2}

    Parameters
    ----------
    V : float
        Molar volume of fluid, [m^3/mol]
    dP_dV : float
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]
    Cp : float
        Real fluid heat capacity at constant pressure, [J/mol/K]
    Cv : float
        Real fluid heat capacity at constant volume, [J/mol/K]
    MW : float, optional
        Molecular weight, [g/mol]

    Returns
    -------
    w : float
        Speed of sound for a real gas, m/s or m*kg^0.5/s/mol^0.5 if MW missing

    Notes
    -----
    An alternate expression based on molar density is as follows:

    .. math::
       w = \left[\left(\frac{\partial P}{\partial \rho}\right)_T \frac{C_p}
       {C_v}\right]^{1/2}

    The form with the unit conversion performed inside it is as follows:

    .. math::
        w = \left[-V^2 \frac{1000}{MW}\left(\frac{\partial P}{\partial V}
        \right)_T \frac{C_p}{C_v}\right]^{1/2}

    Examples
    --------
    Example from [2]_:

    >>> speed_of_sound(V=0.00229754, dP_dV=-3.5459e+08, Cp=153.235, Cv=132.435, MW=67.152)
    179.5868138460819

    References
    ----------
    .. [1] Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey.
       Chemical Thermodynamics for Process Simulation. 1st edition. Weinheim:
       Wiley-VCH, 2012.
    .. [2] Pratt, R. M. "Thermodynamic Properties Involving Derivatives: Using
       the Peng-Robinson Equation of State." Chemical Engineering Education 35,
       no. 2 (March 1, 2001): 112-115.
    '''
    if MW is None:
        return (-V*V*dP_dV*Cp/Cv)**0.5
    else:
        return (-V*V*1000.0*dP_dV*Cp/(Cv*MW))**0.5

def molar_velocity_to_velocity(v_molar, MW):
    r'''Calculate the mass-based velocity (m/s) from the molar velocity of
    the fluid.

    .. math::
        v = \frac{v_{molar}\sqrt{1000}}{\sqrt{\text{MW}}}

    Parameters
    ----------
    v_molar : float
        Molar velcoity, [m*kg^0.5/s/mol^0.5]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    v : float
        Velocity, [m/s]

    Examples
    --------
    >>> molar_velocity_to_velocity(46., 40.445)
    228.73

    '''
    return sqrt(1000)*v_molar/sqrt(MW)

def velocity_to_molar_velocity(v, MW):
    r'''Calculate the molar velocity from the mass-based (m/s) velocity of
    the fluid.

    .. math::
        v_{molar} = \frac{v \sqrt{\text{MW}}}{\sqrt{1000}}

    Parameters
    ----------
    v : float
        Velocity, [m/s]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    v_molar : float
        Molar velcoity, [m*kg^0.5/s/mol^0.5]

    Examples
    --------
    >>> velocity_to_molar_velocity(228.73, 40.445)
    46.
    '''
    return v*sqrt(MW)/sqrt(1000)

def Joule_Thomson(T, V, Cp, dV_dT=None, beta=None):
    r'''Calculate a real fluid's Joule Thomson coefficient. The required
    derivative should be calculated with an equation of state, and `Cp` is the
    real fluid versions. This can either be calculated with `dV_dT` directly,
    or with `beta` if it is already known.

    .. math::
        \mu_{JT} = \left(\frac{\partial T}{\partial P}\right)_H = \frac{1}{C_p}
        \left[T \left(\frac{\partial V}{\partial T}\right)_P - V\right]
        = \frac{V}{C_p}\left(\beta T-1\right)

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    V : float
        Molar volume of fluid, [m^3/mol]
    Cp : float
        Real fluid heat capacity at constant pressure, [J/mol/K]
    dV_dT : float, optional
        Derivative of `V` with respect to `T`, [m^3/mol/K]
    beta : float, optional
        Isobaric coefficient of a thermal expansion, [1/K]

    Returns
    -------
    mu_JT : float
        Joule-Thomson coefficient [K/Pa]

    Examples
    --------
    Example from [2]_:

    >>> Joule_Thomson(T=390, V=0.00229754, Cp=153.235, dV_dT=1.226396e-05)
    1.621956080529905e-05

    References
    ----------
    .. [1] Walas, Stanley M. Phase Equilibria in Chemical Engineering.
       Butterworth-Heinemann, 1985.
    .. [2] Pratt, R. M. "Thermodynamic Properties Involving Derivatives: Using
       the Peng-Robinson Equation of State." Chemical Engineering Education 35,
       no. 2 (March 1, 2001): 112-115.
    '''
    if dV_dT is not None:
        return (T*dV_dT - V)/Cp
    elif beta is not None:
        return V/Cp*(beta*T - 1.)
    else:
        raise ValueError('Either dV_dT or beta is needed')


def isentropic_exponent(Cp, Cv):
    r'''Calculate the isentropic coefficient of an ideal gas, given its constant-
    pressure and constant-volume heat capacity.

    .. math::
        k = \frac{C_p}{C_v}

    Parameters
    ----------
    Cp : float
        Ideal gas heat capacity at constant pressure, [J/mol/K]
    Cv : float
        Ideal gas heat capacity at constant volume, [J/mol/K]

    Returns
    -------
    k : float
        Isentropic exponent, [-]

    Examples
    --------
    >>> isentropic_exponent(33.6, 25.27)
    1.329639889196676

    Notes
    -----
    For real gases, there are more complexities and formulas. Each of the
    formulas reverts to this formula in the case of an ideal gas.

    See Also
    --------
    isentropic_exponent_PV
    isentropic_exponent_PT
    isentropic_exponent_TV

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return Cp/Cv

def isentropic_exponent_PV(Cp, Cv, Vm, P, dP_dV_T):
    r'''Calculate the isentropic coefficient of real fluid using the definition
    of :math:`PV^k = \text{const}`.

    .. math::
        k = -\frac{V}{P}\frac{C_p}{C_v}\left(\frac{\partial P}{\partial V}\right)_T

    Parameters
    ----------
    Cp : float
        Real heat capacity at constant pressure, [J/mol/K]
    Cv : float
        Real heat capacity at constant volume, [J/mol/K]
    Vm : float
        Molar volume, [m^3/mol]
    P : float
        Pressure [Pa]
    dP_dV_T : float
        Derivative of `P` with respect to `V` (at constant temperature),
        [Pa*mol/m^3]

    Returns
    -------
    k_PV : float
        Isentropic exponent of a real fluid, [-]

    Examples
    --------
    Isentropic exponent of air according to Lemmon (2000) at 1000 bar and 300 K:

    >>> isentropic_exponent_PV(Cp=38.36583283578205, Cv=23.98081290153672, Vm=4.730885141495376e-05, P=100000000.0, dP_dV_T=-5417785576072.434)
    4.100576762582646

    See Also
    --------
    isentropic_exponent
    isentropic_exponent_PT
    isentropic_exponent_TV

    References
    ----------
    .. [1] Pini, Matteo. "NiceProp: An Interactive Python-Based Educational
       Tool for Non-Ideal Compressible Fluid Dynamics." SoftwareX 17 (2022):
       100897.
    .. [2] Kouremenos, D. A., and K. A. Antonopoulos. "Isentropic Exponents of
       Real Gases and Application for the Air at Temperatures from 150 K to 450
       K." Acta Mechanica 65, no. 1 (January 1, 1987): 81-99.
       https://doi.org/10.1007/BF01176874.
    '''
    return -Vm*Cp*dP_dV_T/(P*Cv)

def isentropic_exponent_PT(Cp, P, dV_dT_P):
    r'''Calculate the isentropic coefficient of real fluid using the definition
    of :math:`P^{(1-k)}T^k = \text{const}`.

    .. math::
        k = \frac{1}{1 - \frac{P}{C_p}\left(\frac{\partial V}{\partial T}\right)_P}

    Parameters
    ----------
    Cp : float
        Real heat capacity at constant pressure, [J/mol/K]
    P : float
        Pressure [Pa]
    dV_dT_P : float
        Derivative of `V` with respect to `T` (at constant pressure),
        [m^3/(mol*K)]

    Returns
    -------
    k_PT : float
        Isentropic exponent of a real fluid, [-]

    Examples
    --------
    Isentropic exponent of air according to Lemmon (2000) at 1000 bar and 300 K:

    >>> isentropic_exponent_PT(Cp=38.36583283578205, P=100000000.0, dV_dT_P=9.407705210161724e-08)
    1.32487270350443

    See Also
    --------
    isentropic_exponent_PV
    isentropic_exponent
    isentropic_exponent_TV

    References
    ----------
    .. [1] Pini, Matteo. "NiceProp: An Interactive Python-Based Educational
       Tool for Non-Ideal Compressible Fluid Dynamics." SoftwareX 17 (2022):
       100897.
    .. [2] Kouremenos, D. A., and K. A. Antonopoulos. "Isentropic Exponents of
       Real Gases and Application for the Air at Temperatures from 150 K to 450
       K." Acta Mechanica 65, no. 1 (January 1, 1987): 81-99.
       https://doi.org/10.1007/BF01176874.
    '''
    return -Cp/(P*dV_dT_P - Cp) # Avoids a division
    # return 1.0/(1.0 - P*dV_dT_P/Cp)


def isentropic_exponent_TV(Cv, Vm, dP_dT_V):
    r'''Calculate the isentropic coefficient of real fluid using the definition
    of :math:`TV^{k-1} = \text{const}`.

    .. math::
        k = 1 + \frac{V}{C_v} \left(\frac{\partial P}{\partial T}\right)_V

    Parameters
    ----------
    Cv : float
        Real heat capacity at constant volume, [J/mol/K]
    Vm : float
        Molar volume, [m^3/mol]
    dP_dT_V : float
        Derivative of `P` with respect to `T` (at constant volume),
        [Pa/K]

    Returns
    -------
    k_TV : float
        Isentropic exponent of a real fluid, [-]

    Examples
    --------
    Isentropic exponent of air according to Lemmon (2000) at 1000 bar and 300 K:

    >>> isentropic_exponent_TV(Cv=23.98081290153672, Vm=4.730885141495376e-05, dP_dT_V=509689.2959155567)
    2.005504495083

    See Also
    --------
    isentropic_exponent_PV
    isentropic_exponent_PT
    isentropic_exponent

    References
    ----------
    .. [1] Pini, Matteo. "NiceProp: An Interactive Python-Based Educational
       Tool for Non-Ideal Compressible Fluid Dynamics." SoftwareX 17 (2022):
       100897.
    .. [2] Kouremenos, D. A., and K. A. Antonopoulos. "Isentropic Exponents of
       Real Gases and Application for the Air at Temperatures from 150 K to 450
       K." Acta Mechanica 65, no. 1 (January 1, 1987): 81-99.
       https://doi.org/10.1007/BF01176874.
    '''
    return 1.0 + Vm*dP_dT_V/Cv

def Vm_to_rho(Vm, MW):
    r'''Calculate the density of a chemical, given its molar volume and
    molecular weight.

    .. math::
        \rho = \frac{MW}{1000\cdot VM}

    Parameters
    ----------
    Vm : float
        Molar volume, [m^3/mol]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    rho : float
        Density, [kg/m^3]

    Examples
    --------
    >>> Vm_to_rho(0.000132, 86.18)
    652.8787878787879

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return 1e-3*MW/Vm


def rho_to_Vm(rho, MW):
    r'''Calculate the molar volume of a chemical, given its density and
    molecular weight.

    .. math::
        V_m = \left(\frac{1000 \rho}{MW}\right)^{-1}

    Parameters
    ----------
    rho : float
        Density, [kg/m^3]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    Vm : float
        Molar volume, [m^3/mol]

    Examples
    --------
    >>> rho_to_Vm(652.9, 86.18)
    0.0001319957114412621

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return 1e-3*MW/rho


def Z(T, P, V):
    r'''Calculates the compressibility factor of a gas, given its
    temperature, pressure, and molar volume.

    .. math::
        Z = \frac{PV}{RT}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure [Pa]
    V : float
        Molar volume, [m^3/mol]

    Returns
    -------
    Z : float
        Compressibility factor, [-]

    Examples
    --------
    >>> Z(600, P=1E6, V=0.00463)
    0.9281016730797026

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return V*P/(R*T)

def zs_to_ws(zs, MWs):
    r'''Converts a list of mole fractions to mass fractions. Requires molecular
    weights for all species.

    .. math::
        w_i = \frac{z_i MW_i}{MW_{avg}}

    .. math::
        MW_{avg} = \sum_i z_i MW_i

    Parameters
    ----------
    zs : iterable
        Mole fractions [-]
    MWs : iterable
        Molecular weights [g/mol]

    Returns
    -------
    ws : iterable
        Mass fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Examples
    --------
    >>> zs_to_ws([0.5, 0.5], [10, 20])
    [0.3333333333333333, 0.6666666666666666]
    '''
    # Cannot use sum and list comprehension with numba; otherwise Mavg = 1.0/sum(ws)
    # use [0.0]*N initialization for easy transformation into a numpy array in numba
    N = len(zs)
    Mavg = 0.0
    ws = [0.0]*N
    for i in range(N):
        v = zs[i]*MWs[i]
        ws[i] = v
        Mavg += v
    Mavg = 1.0/Mavg
    for i in range(N):
        ws[i] *= Mavg
    return ws

def ws_to_zs(ws, MWs):
    r'''Converts a list of mass fractions to mole fractions. Requires molecular
    weights for all species.

    .. math::
        z_i = \frac{\frac{w_i}{MW_i}}{\sum_i \frac{w_i}{MW_i}}

    Parameters
    ----------
    ws : iterable
        Mass fractions [-]
    MWs : iterable
        Molecular weights [g/mol]

    Returns
    -------
    zs : iterable
        Mole fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Examples
    --------
    >>> ws_to_zs([0.3333333333333333, 0.6666666666666666], [10, 20])
    [0.5, 0.5]
    '''
    N = len(ws)
    zs = [0.0]*N
    tot = 0.0
    for i in range(N):
        v = ws[i]/MWs[i]
        tot += v
        zs[i] = v
    tot = 1.0/tot
    for i in range(N):
        zs[i] *= tot
    return zs


def zs_to_Vfs(zs, Vms):
    r'''Converts a list of mole fractions to volume fractions. Requires molar
    volumes for all species.

    .. math::
        \text{Vf}_i = \frac{z_i V_{m,i}}{\sum_i z_i V_{m,i}}

    Parameters
    ----------
    zs : iterable
        Mole fractions [-]
    Vms : iterable
        Molar volumes of species [m^3/mol]

    Returns
    -------
    Vfs : list
        Molar volume fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Molar volumes are specified in terms of pure components only. Function
    works with any phase.

    Examples
    --------
    Acetone and benzene example

    >>> zs_to_Vfs([0.637, 0.363], [8.0234e-05, 9.543e-05])
    [0.5960229712956298, 0.4039770287043703]
    '''
    N = len(zs)
    Vfs = [0.0]*N
    tot = 0.0
    for i in range(N):
        v = zs[i]*Vms[i]
        tot += v
        Vfs[i] = v
    tot = 1.0/tot
    for i in range(N):
        Vfs[i] *= tot
    return Vfs


def Vfs_to_zs(Vfs, Vms):
    r'''Converts a list of mass fractions to mole fractions. Requires molecular
    weights for all species.

    .. math::
        z_i = \frac{\frac{\text{Vf}_i}{V_{m,i}}}{\sum_i
        \frac{\text{Vf}_i}{V_{m,i}}}

    Parameters
    ----------
    Vfs : iterable
        Molar volume fractions [-]
    Vms : iterable
        Molar volumes of species [m^3/mol]

    Returns
    -------
    zs : list
        Mole fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Molar volumes are specified in terms of pure components only. Function
    works with any phase.

    Examples
    --------
    Acetone and benzene example

    >>> Vfs_to_zs([0.596, 0.404], [8.0234e-05, 9.543e-05])
    [0.6369779395901142, 0.3630220604098858]
    '''
    N = len(Vfs)
    zs = [0.0]*N
    tot = 0.0
    for i in range(N):
        v = Vfs[i]/Vms[i]
        zs[i] = v
        tot += v
    tot = 1.0/tot
    for i in range(N):
        zs[i] *= tot
    return zs

def ms_to_ns(ms, MWs):
    r'''Converts a list of mass flow rates to mole flow rates. Requires molecular
    weights for all species.

    .. math::
        n_i = \frac{1000 m_i}{MW_i}

    Parameters
    ----------
    ms : iterable
        Mass flow rates [kg/s]
    MWs : iterable
        Molecular weights [g/mol]

    Returns
    -------
    ns : iterable
        Mole flow rates [mol/s]

    Notes
    -----
    Does not check that inputs are of the same length.

    Examples
    --------
    >>> ms_to_ns([4, 5], [24, 45])
    [166.666, 111.111]
    '''
    N = len(ms)
    ns = [0.0]*N
    for i in range(N):
        ns[i] = 1e3*ms[i]/MWs[i]
    return ns

def ns_to_ms(ns, MWs):
    r'''Converts a list of mole flow rates to mass flow rates. Requires molecular
    weights for all species.

    .. math::
        m_i = \frac{n_i MW_i}{1000}

    Parameters
    ----------
    ns : iterable
        Mole flow rates [mol/s]
    MWs : iterable
        Molecular weights [g/mol]

    Returns
    -------
    ms : iterable
        Mass flow rates [kg/s]

    Notes
    -----
    Does not check that inputs are of the same length.

    Examples
    --------
    >>> ns_to_ms([166.6666666666, 111.1111111111], [24, 45])
    [4.0, 5.0]
    '''
    N = len(ns)
    ms = [0.0]*N
    for i in range(N):
        ms[i] = ns[i]*MWs[i]*1e-3
    return ms

def ns_to_Qls(ns, Vmls):
    r'''Converts a list of mole flow rates to standard liquid volume
    flow rates. Requires standard liquid molar volumes for all species.

    .. math::
        {Ql}_i = n_i {Vml}_i

    Parameters
    ----------
    ns : iterable
        Mole flow rates [mol/s]
    Vmls : iterable
        Standard molar liquid volumes of each component [m^3/mol]

    Returns
    -------
    Qls : iterable
        Standard liquid volume flow rates [m^3/s]

    Notes
    -----
    Does not check that inputs are of the same length.

    Examples
    --------
    >>> ns_to_Qls([2.0, 3.0], [1e-4, 2e-4])
    [2e-4, 6e-4]
    '''
    N = len(ns)
    Qls = [0.0]*N
    for i in range(N):
        Qls[i] = ns[i] * Vmls[i]
    return Qls

def Qls_to_ns(Qls, Vmls):
    r'''Converts a list of standard liquid volume flow rates
    to mole flow rates. Requires standard liquid molar volumes for all species.

    .. math::
        n_i = \frac{{Ql}_i}{{Vml}_i}

    Parameters
    ----------
    Qls : iterable
        Standard liquid volume flow rates [m^3/s]
    Vmls : iterable
        Standard liquid molar volumes of each component [m^3/mol]

    Returns
    -------
    ns : iterable
        Mole flow rates [mol/s]

    Notes
    -----
    Does not check that inputs are of the same length.

    Examples
    --------
    >>> Qls_to_ns([2e-4, 6e-4], [1e-4, 2e-4])
    [2.0, 3.0]
    '''
    N = len(Qls)
    ns = [0.0]*N
    for i in range(N):
        ns[i] = Qls[i]/Vmls[i]
    return ns

def ms_to_Qls(ms, MWs, Vmls):
    r'''Converts a list of mass flow rates to standard liquid volume
    flow rates. Requires molecular weights and standard molar liquid
    volumes for all species.

    .. math::
        {Ql}_i = \frac{1000 m_i {Vml}_i}{MW_i}

    Parameters
    ----------
    ms : iterable
        Mass flow rates [kg/s]
    MWs : iterable
        Molecular weights [g/mol]
    Vmls : iterable
        Standard liquid molar volumes [m^3/mol]

    Returns
    -------
    Qls : iterable
        Standard liquid volume flow rates [m^3/s]

    Notes
    -----
    Does not check that inputs are of the same length.

    Examples
    --------
    >>> ms_to_Qls([4.0, 5.0], [24, 45], [1e-4, 2e-4])
    [0.0166666666, 0.0222222222]
    '''
    N = len(ms)
    Qls = [0.0]*N
    for i in range(N):
        Qls[i] = 1e3*ms[i]/MWs[i]*Vmls[i]
    return Qls

def Qls_to_ms(Qls, MWs, Vmls):
    r'''Converts a list of standard liquid volume flow rates to mass
    flow rates. Requires molecular weights and standard liquid molar
    volumes for all species.

    .. math::
        m_i = \frac{{Ql}_i {MW}_i}{1000 {Vml}_i}

    Parameters
    ----------
    Qls : iterable
        Standard liquid volume flow rates [m^3/s]
    MWs : iterable
        Molecular weights [g/mol]
    Vmls : iterable
        Molar volumes in the liquid phase [m^3/mol]

    Returns
    -------
    ms : iterable
        Mass flow rates [kg/s]

    Notes
    -----
    Does not check that inputs are of the same length.

    Examples
    --------
    >>> Qls_to_ms([1.666666666e-02, 1.11111111e-01], [24, 45], [1e-4, 2e-4])
    [4.0, 25.0]
    '''
    N = len(Qls)
    ms = [0.0]*N
    for i in range(N):
        ms[i] = 1e-3*Qls[i]*MWs[i]/Vmls[i]
    return ms

def dxs_to_dns(dxs, xs, dns=None):
    r'''Convert the mole fraction derivatives of a quantity (calculated so
    they do not sum to 1) to mole number derivatives (where the mole fractions
    do sum to one). Requires the derivatives and the mole fractions of the
    mixture.

    .. math::
       \left(\frac{\partial M}{\partial n_i}\right)_{n_{k\ne i}} = \left[
        \left(\frac{\partial M}{\partial x_i}\right)_{x_{k\ne i}}
        - \sum_j x_j \left(\frac{\partial M}{\partial x_j}  \right)_{x_{k\ne j}}
        \right]

    Parameters
    ----------
    dxs : list[float]
        Derivatives of a quantity with respect to mole fraction (not summing to
        1), [prop]
    xs : list[float]
        Mole fractions of the species, [-]
    dns : list[float], optional
        Return array, [prop/mol]

    Returns
    -------
    dns : list[float]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop/mol]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    This applies to a specific phase only, not to a mixture of multiple phases.

    Examples
    --------
    >>> dxs_to_dns([-0.0028, -0.00719, -0.00859], [0.7, 0.2, 0.1])
    [0.0014570000000000004, -0.002933, -0.004333]
    '''
    xdx_tot = 0.0
    N = len(xs)
    for j in range(N):
        xdx_tot += xs[j]*dxs[j]
    if dns is None:
        dns = [0.0]*N
    for i in range(N):
        dns[i] = dxs[i] - xdx_tot
    return dns


def dns_to_dn_partials(dns, F, partial_properties=None):
    r'''Convert the mole number derivatives of a quantity (calculated so
    they do sum to 1) to partial molar quantites.

    .. math::

        \left(\frac{\partial n F}{\partial n_i}\right)_{n_{k \ne i}} = F_i +
        n \left(\frac{\partial F}{\partial n_i}\right)_{n_{k\ne i}}

    In the formula, the `n` is 1.

    Parameters
    ----------
    dns : list[float]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop/mol]
    F : float
        Property evaluated at constant composition, [prop]
    partial_properties : list[float], optional
        Optional output array for derivatives of a quantity with respect
        to mole number (summing to 1), [prop]

    Returns
    -------
    partial_properties : list[float]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    This applies to a specific phase only, not to a mixture of multiple phases.

    This is especially useful for fugacity calculations.

    Examples
    --------
    >>> dns_to_dn_partials([0.001459, -0.002939, -0.004334], -0.0016567)
    [-0.0001977000000000001, -0.0045957, -0.0059907]
    '''
    N = len(dns)
    if partial_properties is None:
        partial_properties = [0.0]*N
    for i in range(N):
        partial_properties[i] = F + dns[i]
    return partial_properties


def dxs_to_dn_partials(dxs, xs, F, partial_properties=None):
    r'''Convert the mole fraction derivatives of a quantity (calculated so
    they do not sum to 1) to partial molar quantites. Requires the derivatives
    and the mole fractions of the mixture.

    .. math::

        \left(\frac{\partial n F}{\partial n_i}\right) =
        \left(\frac{\partial F}{\partial x_i}\right)+ F
        - \sum_j x_j \left(\frac{\partial F}{\partial x_j}\right)


    Parameters
    ----------
    dxs : list[float]
        Derivatives of a quantity with respect to mole fraction (not summing to
        1), [prop]
    xs : list[float]
        Mole fractions of the species, [-]
    F : float
        Property evaluated at constant composition, [prop]
    partial_properties : list[float], optional
        Array for Derivatives of a quantity with respect to mole number (summing to
        1), [prop]

    Returns
    -------
    partial_properties : list[float]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    This applies to a specific phase only, not to a mixture of multiple phases.

    See Also
    --------
    dxs_to_dns
    dns_to_dn_partials

    Examples
    --------
    >>> dxs_to_dn_partials([-0.0026404, -0.00719, -0.00859], [0.7, 0.2, 0.1],
    ... -0.0016567)
    [-0.00015182, -0.0047014199999999996, -0.00610142]
    '''
    xdx_totF = F
    N = len(xs)
    if partial_properties is None:
        partial_properties = [0.0]*N
    for j in range(N):
        xdx_totF -= xs[j]*dxs[j]
    for i in range(N):
        partial_properties[i] = dxs[i] + xdx_totF
    return partial_properties


def d2ns_to_dn2_partials(d2ns, dns):
    r'''Convert second-order mole number derivatives of a quantity
     to the following second-order partial derivative:

    .. math::
            \frac{\partial^2 n F}{\partial n_j \partial n_i}
            = \frac{\partial^2 F}{\partial n_i \partial n_j}
            + \frac{\partial F}{\partial n_i}
            + \frac{\partial F}{\partial n_j}


    Requires the second order mole number derivatives and the first order
    mole number derivatives of the mixture only.

    Parameters
    ----------
    d2ns : list[float]
        Second order derivatives of a quantity with respect to mole number
        (summing to 1), [prop/mol^2]
    dns : list[float]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop/mol]

    Returns
    -------
    second_partial_properties : list[list[float]]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    This was originally implemented to allow for the calculation of
    first mole number derivatices of log fugacity coefficients; the two
    arguments are the second and first mole number derivatives of the overall
    mixture log fugacity coefficient.


    Derived with the following SymPy code.

    >>> from sympy import * # doctest: +SKIP
    >>> n1, n2 = symbols('n1, n2') # doctest: +SKIP
    >>> f, g, h = symbols('f, g, h', cls=Function) # doctest: +SKIP
    >>> diff(h(n1, n2)*f(n1,  n2), n1, n2) # doctest: +SKIP
    f(n1, n2)*Derivative(h(n1, n2), n1, n2) + h(n1, n2)*Derivative(f(n1, n2), n1, n2) + Derivative(f(n1, n2), n1)*Derivative(h(n1, n2), n2) + Derivative(f(n1, n2), n2)*Derivative(h(n1, n2), n1)

    See Also
    --------
    dxs_to_dns
    dns_to_dn_partials
    dxs_to_dn_partials

    Examples
    --------
    >>> d2ns = [[0.152, 0.08, 0.547], [0.08, 0.674, 0.729], [0.547, 0.729, 0.131]]
    >>> d2ns_to_dn2_partials(d2ns, [20.0, .124, 900.52])
    [[40.152, 20.203999999999997, 921.067], [20.204, 0.922, 901.3729999999999], [921.067, 901.373, 1801.1709999999998]]
    '''
    cmps = range(len(dns))
    hess = []
    for i in cmps:
        row = []
        for j in cmps:
            v = d2ns[i][j] + dns[i] + dns[j]
            row.append(v)
        hess.append(row)
    return hess



def d2xs_to_dxdn_partials(d2xs, xs):
    r'''Convert second-order mole fraction derivatives of a quantity
    (calculated so they do not sum to 1) to the following second-order
    partial derivative:

    .. math::
            \frac{\partial^2 n F}{\partial x_j \partial n_i}
            = \frac{\partial^2 F}{\partial x_i x_j}
            - \sum_k x_k \frac{\partial^2 F}{\partial x_k \partial x_j}

    Requires the second derivatives and the mole fractions of the mixture only.

    Parameters
    ----------
    d2xs : list[float]
        Second derivatives of a quantity with respect to mole fraction (not
         summing to 1), [prop]
    xs : list[float]
        Mole fractions of the species, [-]

    Returns
    -------
    partial_properties : list[float]
        Derivatives of a quantity with respect to mole number (summing to
        1), [prop]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    See Also
    --------
    dxs_to_dns
    dns_to_dn_partials
    dxs_to_dn_partials

    Examples
    --------
    >>> d2xs = [[0.152, 0.08, 0.547], [0.08, 0.674, 0.729], [0.547, 0.729, 0.131]]
    >>> d2xs_to_dxdn_partials(d2xs, [0.7, 0.2, 0.1])
    [[-0.02510000000000001, -0.18369999999999997, 0.005199999999999982], [-0.0971, 0.41030000000000005, 0.18719999999999992], [0.3699, 0.4653, -0.41080000000000005]]
    '''
    N = len(xs)

    double_sums = [0.0]*N
    for j in range(N):
        tot = 0.0
        for k in range(N):
            tot += xs[k]*d2xs[j][k]
        double_sums[j] = tot

    # Oddly, the below saw found to be successful for NRTL but not PR
    # Mysterious, interesting sum which is surprisingly efficient to calculate
#    d2xsi = d2xs[1]
#    symmetric_sum = 0.0
#    for k in cmps:
#        symmetric_sum += xs[k]*d2xsi[k]
#    print(symmetric_sum)

    return [[d2xs[i][j] - double_sums[j] for j in range(N)]
             for i in range(N)]



def dxs_to_dxsn1(dxs):
    r'''Convert the mole fraction derivatives of a quantity (calculated so
    they do not sum to 1) to derivatives such that they do sum to 1 by changing
    the composition of the last component in the negative of the component
    which is changed. Requires the derivatives of the mixture only. The size of
    the returned array is one less than the input.

    .. math::
        \left(\frac{\partial F}{\partial x_i}\right)_{\sum_{x_i}^N =1} =
        \left(\frac{\partial F}{\partial x_i}
        - \frac{\partial F}{\partial x_N}\right)_{\sum_{x_i}^N \ne 1}

    Parameters
    ----------
    dxs : list[float]
        Derivatives of a quantity with respect to mole fraction (not summing to
        1), [prop]

    Returns
    -------
    dxsm1 : list[float]
        Derivatives of a quantity with respect to mole fraction (summing to
        1 by altering the last component's composition), [prop]

    Notes
    -----

    Examples
    --------
    >>> dxs_to_dxsn1([-2651.3181821109024, -2085.574403592012, -2295.0860830203587])
    [-356.23209909054367, 209.51167942834672]
    '''
    last = dxs[-1]
    return [dx - last for dx in dxs[:-1]]


def d2xs_to_d2xsn1(d2xs):
    r'''Convert the second mole fraction derivatives of a quantity (calculated
    so they do not sum to 1) to derivatives such that they do sum to 1
    Requires the second derivatives of the mixture only. The size of
    the returned array is one less than the input in both dimensions

    .. math::
        \left(\frac{\partial^2 F}{\partial x_i \partial x_j }\right)_{\sum_{x_i}^N =1} =
        \left(\frac{\partial^2 F}{\partial x_i\partial x_j}
        -\frac{\partial^2 F}{\partial x_i\partial x_N}
        -\frac{\partial^2 F}{\partial x_j\partial x_N}
        +\frac{\partial^2 F}{\partial x_N\partial x_N}
        \right)_{\sum_{x_i}^N \ne 1}

    Parameters
    ----------
    d2xs : list[float]
        Second of a quantity with respect to mole fraction (not summing to
        1), [prop]

    Returns
    -------
    d2xsm1 : list[float]
        Second derivatives of a quantity with respect to mole fraction (summing
        to 1 by altering the last component's composition), [prop]

    Notes
    -----

    Examples
    --------
    >>> d2xs_to_d2xsn1([[-2890.4327598108, -6687.0990540960065, -1549.375443699441], [-6687.099054095983, -2811.2832904869883, -1228.6223853777503], [-1549.3754436994498, -1228.6223853777562, -3667.388098758508]])
    [[-3459.069971170426, -7576.489323777324], [-7576.489323777299, -4021.4266184899957]]
    '''
    N = len(d2xs)-1
    out = [[0.0]*N for _ in range(N)]
    last = d2xs[-1][-1]
    for i in range(N):
        last_i = d2xs[i][-1]
        for j in range(N):
            out[i][j] = d2xs[i][j] - last_i - d2xs[j][-1] + last
    return out


def none_and_length_check(all_inputs, length=None):
    r'''Checks inputs for suitability of use by a mixing rule which requires
    all inputs to be of the same length and non-None. A number of variations
    were attempted for this function; this was found to be the quickest.

    Parameters
    ----------
    all_inputs : array-like of array-like
        list of all the lists of inputs, [-]
    length : int, optional
        Length of the desired inputs, [-]

    Returns
    -------
    False/True : bool
        Returns True only if all inputs are the same length (or length `length`)
        and none of the inputs contain None [-]

    Notes
    -----
    Does not check for nan values.

    Examples
    --------
    >>> none_and_length_check(([1, 1], [1, 1], [1, 30], [10,0]), length=2)
    True
    '''
    if not length:
        length = len(all_inputs[0])
    return all(not (None in things or len(things) != length) for things in all_inputs)




def normalize(values):
    r'''Simple function which normalizes a series of values to be from 0 to 1,
    and for their sum to add to 1.

    .. math::
        x = \frac{x}{sum_i x_i}

    Parameters
    ----------
    values : array-like
        array of values

    Returns
    -------
    fractions : array-like
        Array of values from 0 to 1

    Notes
    -----
    Does not work on negative values, or handle the case where the sum is zero.

    Examples
    --------
    >>> normalize([3, 2, 1])
    [0.5, 0.3333333333333333, 0.16666666666666666]
    '''
    try:
        tot_inv = 1.0/sum(values)
        return [i*tot_inv for i in values]
    except ZeroDivisionError:
        N = len(values)
        try:
            # case of values sum to zero
            return [1.0/N]*N
        except ZeroDivisionError:
            # case of 0 values
            return []

def remove_zeros(values, tol=1e-6):
    r'''Simple function which removes zero values from an array, and replaces
    them with a user-specified value, normally a very small number. Helpful
    for the case where a function can work with values very close to zero but
    not quite zero. The resulting array is normalized so the sum is still one.

    Parameters
    ----------
    values : array-like
        array of values
    tol : float
        The replacement value for zeroes

    Returns
    -------
    fractions :  array-like
        Array of values from 0 to 1

    Notes
    -----
    Works on numpy arrays, and returns numpy arrays only for that case.

    Examples
    --------
    >>> remove_zeros([0, 1e-9, 1], 1e-12)
    [9.99999998999e-13, 9.99999998999e-10, 0.999999998999]
    '''
    if any(i == 0 for i in values):
        ans = normalize([i if i != 0 else tol for i in values])
        if type(values) != list:
            if isinstance(values, np.ndarray):
                return np.array(ans)
        return ans
    return values


def mixing_simple(fracs, props):
    r'''Simple function calculates a property based on weighted averages of
    properties. Weights could be mole fractions, volume fractions, mass
    fractions, or anything else.

    .. math::
        y = \sum_i \text{frac}_i \cdot \text{prop}_i

    Parameters
    ----------
    fracs : array-like
        Fractions of a mixture
    props : array-like
        Properties

    Returns
    -------
    prop : value
        Calculated property

    Notes
    -----

    Examples
    --------
    >>> mixing_simple([0.1, 0.9], [0.01, 0.02])
    0.019000000000000003
    '''
    tot = 0.0
    for i in range(len(fracs)):
        tot += fracs[i]*props[i]
    return tot


def mixing_logarithmic(fracs, props):
    r'''Simple function calculates a property based on weighted averages of
    logarithmic properties.

    .. math::
        y = \sum_i \text{frac}_i \cdot \ln(\text{prop}_i)

    Parameters
    ----------
    fracs : array-like
        Fractions of a mixture
    props : array-like
        Properties

    Returns
    -------
    prop : value
        Calculated property

    Notes
    -----
    Does not work on negative values.
    Returns None if any fractions or properties are missing or are not of the
    same length.

    Examples
    --------
    >>> mixing_logarithmic([0.1, 0.9], [0.01, 0.02])
    0.01866065983073615
    '''
    tot = 0.0
    for i in range(len(fracs)):
        tot += fracs[i]*trunc_log(props[i])
    return trunc_exp(tot)

def mixing_power(fracs, props, r):
    r'''Power law mixing rule for any property, with a variable exponent
    `r` as input. Optimiezd routines are available for common powers.

    .. math::
        \text{prop}_{mix}^r = \sum_i z_i \left(\text{prop}_i \right)^{r}

    Parameters
    ----------
    fracs : list[float]
        Mole fractions of components (or mass, or volume, etc.), [-]
    props : list[float]
        Properties of all components, [various]
    r : float
        Power mixing exponent, [-]

    Returns
    -------
    prop : float
        Property for mixture, [`props`]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.


    The following recommendations in [1] exist for different properties:

    **Surface tension**: r = 1 Recommended by an author in [1]_; but often
    non-linear behavior is shown and r= -1 to r=-3 is recommended.  r = -1
    is most often used.

    **Liquid thermal conductivity**: r = -2 in [1]_; this is known also as
    procedure DIPPR9B.

    Examples
    --------
    >>> mixing_power([0.258, 0.742], [0.1692, 0.1528], -2)
    0.15657104706719646

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    N = len(fracs)
    prop = 0.0
    if r == -4.0:
        for i in range(N):
            prop_i_squared = props[i]*props[i]
            prop += fracs[i]/(prop_i_squared*prop_i_squared)
        return prop**(1.0/r)
    elif r == -3.0:
        for i in range(N):
            prop += fracs[i]/(props[i]*props[i]*props[i])
        return prop**(1.0/r)
    elif r == -2.0:
        for i in range(N):
            prop += fracs[i]/(props[i]*props[i])
        return 1.0/sqrt(prop)
    elif r == -1.0:
        for i in range(N):
            prop += fracs[i]/(props[i])
        return 1.0/(prop)
    elif r == 1.0:
        for i in range(N):
            prop += fracs[i]*(props[i])
        return (prop)
    elif r == 2.0:
        for i in range(N):
            prop += fracs[i]*(props[i]*props[i])
        return sqrt(prop)
    elif r == 3.0:
        for i in range(N):
            prop += fracs[i]*(props[i]*props[i]*props[i])
        return prop**(1.0/3.0)
    elif r == 4.0:
        for i in range(N):
            prop_i_squared = props[i]*props[i]
            prop += fracs[i]*(prop_i_squared*prop_i_squared)
        return sqrt(sqrt(prop))
    elif r == 1/3:
        for i in range(N):
            prop += fracs[i] * cbrt(props[i])
        return prop * prop * prop
    elif r == -1/3:
        for i in range(N):
            prop += fracs[i] / cbrt(props[i])
        return 1.0 / (prop * prop * prop)
    elif r == 2/3:
        for i in range(N):
            prop += fracs[i] * cbrt(props[i] * props[i])
        return prop * sqrt(prop)
    elif r == -2/3:
        for i in range(N):
            prop += fracs[i] / cbrt(props[i] * props[i])
        return 1.0 / (prop * sqrt(prop))
    elif r == 4/3:
        for i in range(N):
            prop += fracs[i] * (props[i] * cbrt(props[i]))
        rt_prop = sqrt(prop)
        return sqrt(rt_prop)*rt_prop
    elif r == -4/3:
        for i in range(N):
            prop += fracs[i]/(props[i] * cbrt(props[i]))
        rt_prop = sqrt(prop)
        return 1.0/(sqrt(rt_prop)*rt_prop)
        
    elif r == 0.5:
        for i in range(N):
            prop += fracs[i] * sqrt(props[i])
        return prop*prop   
    elif r == -0.5:
        for i in range(N):
            prop += fracs[i] / sqrt(props[i])
        return 1.0 / (prop*prop)
    elif r == 1.5:
        for i in range(N):
            prop += fracs[i]*(props[i]*props[i])
        return prop**(2.0/3.0)
    elif r == -1.5:
        for i in range(N):
            prop += fracs[i] / (props[i]*sqrt(props[i]))
        return prop**(-2.0/3.0)
    elif r == 0.25:
        for i in range(N):
            prop += fracs[i] * sqrt(sqrt(props[i]))
        prop2 = prop*prop
        return prop2*prop2
    elif r == -0.25:
        for i in range(N):
            prop += fracs[i] / sqrt(sqrt(props[i]))
        prop2 = prop*prop
        return 1.0 / (prop2*prop2)
    elif r == 0.75:
        for i in range(N):
            rt_prop_i = sqrt(props[i])
            prop += fracs[i] * (sqrt(rt_prop_i)*rt_prop_i)
        return prop ** (4/3)
    elif r == -0.75:
        for i in range(N):
            rt_prop_i = sqrt(props[i])
            prop += fracs[i] / (sqrt(rt_prop_i)*rt_prop_i)
        return (prop ** (-4/3))
    elif r == -0.125:
        for i in range(N):
            prop += fracs[i] / sqrt(sqrt(sqrt(props[i])))
        prop *= prop
        prop *= prop
        return 1.0 / (prop*prop)
    elif r == 0.125:
        for i in range(N):
            prop += fracs[i] * sqrt(sqrt(sqrt(props[i])))
        prop *= prop
        prop *= prop
        return prop*prop
    elif r == 0.875:
        for i in range(N):
            rt_prop = sqrt(props[i])
            qtrt_prop = sqrt(rt_prop)
            eight_rt_prop = sqrt(qtrt_prop)
            prop += fracs[i] * (rt_prop*qtrt_prop*eight_rt_prop)
        return prop ** (8/7)
    elif r == -0.875:
        for i in range(N):
            rt_prop = sqrt(props[i])
            qtrt_prop = sqrt(rt_prop)
            eight_rt_prop = sqrt(qtrt_prop)
            prop += fracs[i] /(rt_prop*qtrt_prop*eight_rt_prop)
        return prop ** (-8/7)
    elif r == 0.375:
        for i in range(N):
            qtrt_prop = sqrt(sqrt(props[i]))
            eight_rt_prop = sqrt(qtrt_prop)
            prop += fracs[i] * (qtrt_prop*eight_rt_prop)
        return prop ** (1.0/0.375)
    elif r == -0.375:
        for i in range(N):
            qtrt_prop = sqrt(sqrt(props[i]))
            eight_rt_prop = sqrt(qtrt_prop)
            prop += fracs[i] /(qtrt_prop*eight_rt_prop)
        return prop ** (-1.0/0.375)
    elif r == 0.625:
        for i in range(N):
            rt_prop = sqrt(props[i])
            qtrt_prop = sqrt(rt_prop)
            eight_rt_prop = sqrt(qtrt_prop)
            prop += fracs[i] * (rt_prop*eight_rt_prop)
        return prop ** (1.0/0.625)
    elif r == -0.625:
        for i in range(N):
            rt_prop = sqrt(props[i])
            qtrt_prop = sqrt(rt_prop)
            eight_rt_prop = sqrt(qtrt_prop)
            prop += fracs[i] /(rt_prop*eight_rt_prop)
        return prop ** (-1.0/0.625)

    elif r == 5.0/3.0:
        for i in range(N):
            cbrt_prop = cbrt(props[i])
            prop += fracs[i] * (props[i]* cbrt_prop*cbrt_prop)
        return prop**(1.0/(5/3.0))
    elif r == -5.0/3.0:
        for i in range(N):
            cbrt_prop = cbrt(props[i])
            prop += fracs[i] /(props[i]* cbrt_prop*cbrt_prop)
        return prop**(-1.0/(5.0/3.0))
    elif r == 1.0/6.0:
        for i in range(N):
            prop += fracs[i] * sqrt(cbrt(props[i]))
        prop *= prop
        return prop*prop*prop
    elif r == 5/6:
        for i in range(N):
            prop_sqrt = sqrt(props[i])
            prop_1_6 = cbrt(prop_sqrt)
            prop += fracs[i] * (prop_sqrt*prop_1_6*prop_1_6)
        return prop ** (6.0/5.0)
    elif r == -1.0/6.0:
        for i in range(N):
            prop += fracs[i] / sqrt(cbrt(props[i]))
        prop *= prop
        return 1.0/(prop*prop*prop)
    elif r == -5/6:
        for i in range(N):
            prop_sqrt = sqrt(props[i])
            prop_1_6 = cbrt(prop_sqrt)
            prop += fracs[i] /(prop_sqrt*prop_1_6*prop_1_6)
        return prop ** (-6.0/5.0)

    for i in range(N):
        prop += fracs[i]*(props[i]**r)
    return prop**(1.0/r)


def mix_component_flows(IDs1, IDs2, flow1, flow2, fractions1, fractions2):
    r'''Mix two flows of potentially different chemicals of given overall flow
    rates and flow fractions to determine the outlet components, flow rates,
    and compositions. The flows do not need to be of the same length.

    Parameters
    ----------
    IDs1 : list[str]
        List of identifiers of the chemical species in flow one, [-]
    IDs2 : list[str]
        List of identifiers of the chemical species in flow two, [-]
    flow1 : float
        Total flow rate of the chemicals in flow one, [mol/s]
    flow2 : float
        Total flow rate of the chemicals in flow two, [mol/s]
    fractions1 : list[float]
        Mole fractions of each chemical in flow one, [-]
    fractions2 : list[float]
        Mole fractions of each chemical in flow two, [-]

    Returns
    -------
    cmps : list[str]
        List of identifiers of the chemical species in the combined flow, [-]
    moles : list[float]
        Flow rates of all chemical species in the combined flow, [mol/s]

    Notes
    -----
    Mass or volume flows and fractions can be used instead of molar ones.

    If the two flows have the same components, the output list will be in the
    same order as the one given; otherwise they are sorted alphabetically.

    Examples
    --------
    >>> mix_component_flows(['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], 1, 1, [0.5, 0.5], [0.5, 0.5])
    (['64-17-5', '67-56-1', '7732-18-5'], [0.5, 0.5, 1.0])
    '''
    if (set(IDs1) == set(IDs2)) and (len(IDs1) == len(IDs2)):
            cmps = IDs1
    else:
        cmps = sorted(list(set(IDs1 + IDs2)))
    mole = flow1 + flow2
    moles = []
    for cmp in cmps:
        moles.append(0)
        if cmp in IDs1:
            ind = IDs1.index(cmp)
            moles[-1] += fractions1[ind]*flow1
        if cmp in IDs2:
            ind = IDs2.index(cmp)
            moles[-1] += fractions2[ind]*flow2
    return cmps, moles


def mix_component_partial_flows(IDs1, IDs2, ns1=None, ns2=None):
    r'''Mix two flows of potentially different chemicals; with the feature that
    the mole flows of either or both streams may be unknown.

    The flows do not need to be of the same length.

    Parameters
    ----------
    IDs1 : list[str]
        List of identifiers of the chemical species in flow one, [-]
    IDs2 : list[str]
        List of identifiers of the chemical species in flow two, [-]
    ns1 : list[float]
        Total flow rate of the chemicals in flow one, [mol/s]
    ns2 : list[float]
        Total flow rate of the chemicals in flow two, [mol/s]

    Returns
    -------
    cmps : list[str]
        List of identifiers of the chemical species in the combined flow, [-]
    moles : list[float]
        Flow rates of all chemical species in the combined flow, [mol/s]

    Notes
    -----
    Mass or volume flows and fractions can be used instead of molar ones.

    If the two flows have the same components, the output list will be in the
    same order as the one given; otherwise they are sorted alphabetically.

    Examples
    --------
    >>> mix_component_partial_flows(['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], [0.5, 0.5], [0.5, 0.5])
    (['64-17-5', '67-56-1', '7732-18-5'], [0.5, 0.5, 1.0])
    >>> mix_component_partial_flows(['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], None, [0.5, 0.5])
    (['64-17-5', '67-56-1', '7732-18-5'], [0.0, 0.5, 0.5])
    >>> mix_component_partial_flows(['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], [0.5, 0.5], None)
    (['64-17-5', '67-56-1', '7732-18-5'], [0.5, 0.0, 0.5])
    >>> mix_component_partial_flows(['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], None, None)
    (['64-17-5', '67-56-1', '7732-18-5'], [0.0, 0.0, 0.0])
    '''
    if (set(IDs1) == set(IDs2)) and (len(IDs1) == len(IDs2)):
        cmps = IDs1
    else:
        cmps = sorted(list(set(IDs1 + IDs2)))

    ns = []
    for ID in cmps:
        ns.append(0.0)
        if ID in IDs1 and ns1 is not None:
            ind = IDs1.index(ID)
            ns[-1] += ns1[ind]

        if ID in IDs2 and ns2 is not None:
            ind = IDs2.index(ID)
            ns[-1] += ns2[ind]
    return cmps, ns


def mix_multiple_component_flows(IDs, flows, fractions):
    r'''Mix multiple flows of potentially different chemicals of given overall
    flow rates and flow fractions to determine the outlet components, flow
    rates,  and compositions. The flows do not need to be of the same length.

    Parameters
    ----------
    IDs : list[list[str]]
        List of lists of identifiers of the chemical species in the flows, [-]
    flows : list[float]
        List of total flow rates of the chemicals in the streams, [mol/s]
    fractions : list[list[float]]
        List of lists of mole fractions of each chemical in each flow, [-]

    Returns
    -------
    cmps : list[str]
        List of identifiers of the chemical species in the combined flow, [-]
    moles : list[float]
        Flow rates of all chemical species in the combined flow, [mol/s]

    Notes
    -----
    Mass or volume flows and fractions can be used instead of molar ones.

    If the every flow have the same components, the output list will be in the
    same order as the one given; otherwise they are sorted alphabetically.

    Examples
    --------
    >>> mix_multiple_component_flows([['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1']],
    ... [1, 1], [[0.5, 0.5], [0.5, 0.5]])
    (['64-17-5', '67-56-1', '7732-18-5'], [0.5, 0.5, 1.0])
    '''
    n_inputs = len(IDs)
    assert n_inputs == len(flows) == len(fractions)
    if n_inputs == 1:
        n = flows[0]
        return IDs[0], [zi*n for zi in fractions[0]]
    else:
        cmps, component_flows = mix_component_flows(IDs[0], IDs[1], flows[0], flows[1], fractions[0], fractions[1])
        if n_inputs == 2:
            return cmps, component_flows
        else:
            flow = sum(component_flows)
            fracs = [i/flow for i in component_flows]
            counter = 2
            while counter != n_inputs:
                cmps, component_flows = mix_component_flows(cmps, IDs[counter], flow, flows[counter], fracs, fractions[counter])
                flow = sum(component_flows)
                fracs = [i/flow for i in component_flows]
                counter +=1
    return cmps, component_flows




def solve_flow_composition_mix(Fs, zs, ws, MWs):
    r'''Solve a stream composition problem where some specs are mole flow rates;
    some are mass fractions; and some are mole fractions. This algorithm
    requires at least one mole flow rate; and for every other component, a
    single spec in mole or mass or a flow rate. It is permissible for no
    components to have mole fractions; or no components to have weight
    fractions; or both.

    Parameters
    ----------
    Fs : list[float]
        List of mole flow rates; None if not specified for a component, [mol/s]
    zs : list[float]
        Mole fractions; None if not specified for a component [-]
    ws : list[float]
        Mass fractions; None if not specified for a component [-]
    MWs : list[float]
        Molecular weights, [g/mol]

    Returns
    -------
    Fs : list[float]
        List of mole flow rates, [mol/s]
    zs : list[float]
        Mole fractions, [-]
    ws : list[float]
        Mass fractions, [-]

    Notes
    -----
    A fast path is used if no weight fractions are provided; the calculation is
    much simpler for that case.

    This algorithm was derived using SymPy, and framed in a form which allows
    for explicit solving.
    This is capable of solving large-scale problems i.e. with 1000 components a
    solve time is 1 ms; with 10000 it is 10 ms.

    Examples
    --------
    >>> Fs = [3600, None, None, None, None]
    >>> zs = [None, .1, .2, None, None]
    >>> ws = [None, None, None, .01, .02]
    >>> MWs = [18.01528, 46.06844, 32.04186, 72.151, 142.286]
    >>> Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    >>> Fs
    [3600, 519.3039148597746, 1038.6078297195493, 17.44015034881175, 17.687253669610733]
    >>> zs
    [0.6932356751002141, 0.1, 0.2, 0.0033583706669188186, 0.003405954232867038]
    >>> ws
    [0.5154077420893426, 0.19012206531421305, 0.26447019259644433, 0.01, 0.02]
    '''
    # Fs needs to have at least one flow in it
    # Either F, z, or w needs to be specified for every flow.
    # MW needs to be specified for every component.
    zs = list(zs)
    ws = list(ws)
    Fs = list(Fs)

    N = len(MWs)
    F_knowns = [F for F in Fs if F is not None]
    MWs_known = [MWs[i] for i in range(N) if Fs[i] is not None]
    F_known = sum(F_knowns)
    F_known_inv = 1.0/F_known
    MW_known = sum([F*MW*F_known_inv for F, MW in zip(F_knowns, MWs_known)])

    if not any(i is not None for i in ws):
        # Only flow rates or mole fractions!
        zs_spec = sum([zi for zi in zs if zi is not None])
        F_act = F_known/(1.0 - zs_spec)
        for i in range(len(Fs)):
            if Fs[i] is None:
                Fs[i] = zs[i]*F_act

        Ft_inv = 1.0/sum(Fs)
        zs = [F*Ft_inv for F in Fs]
        ws = zs_to_ws(zs, MWs)
        return Fs, zs, ws

    den_bracketed = sum([w for w in ws if w is not None]) - 1.0 # Move out of loop

    MW_z_sum = sum([MWs[i]*zs[i] for i in range(N) if zs[i] is not None])
    zs_sum = sum([zs[i] for i in range(N) if zs[i] is not None])
    weight_term = sum([-ws[i]/(MWs[i]*den_bracketed) for i in range(N) if ws[i] is not None])


    Ft = -F_known*(MW_known*weight_term + 1.0)/(MW_z_sum*weight_term + zs_sum - 1.0)
    Ft_inv = 1.0/Ft
    # Once the problem had been reduced by putting all variables outside the
    # main loop, it was easy to get an analytical solution
#    def to_solve(Ft):
#        num_bracketed = MW_z_sum*Ft + F_known*MW_known
#        F_weights = weight_term*num_bracketed
#        F_fracs = zs_sum*Ft
#        err = (F_known + F_fracs + F_weights) - Ft
#        return err
#    Ft = newton(to_solve, F_known)
    # For each flow, in-place replace with mole flows calculated

    num_bracketed = MW_z_sum*Ft + F_known*MW_known

    for i in range(N):
        if zs[i] is not None:
            Fs[i] = Ft*zs[i]
        if ws[i] is not None:
            Fs[i] = -ws[i]*num_bracketed/(MWs[i]*den_bracketed)

    zs = [F*Ft_inv for F in Fs]
    ws = zs_to_ws(zs, MWs)

    return Fs, zs, ws



def radius_of_gyration(MW, A, B, C, planar=False):
    r'''Calculates the radius of gyration of a molecule using the DIPPR
    definition. The parameters `A`, `B`, and `C` must be obtained from
    either vibrational scpectra and analysis or quantum chemistry calculations
    of programs such as `psi <https://psicode.org/>`.

    For planar molecules defined by only two moments of inertia,

    .. math::
        R_g = \sqrt{\sqrt{AB}\frac{N_A}{\text{MW}}}

    For non-planar molecules with three moments of inertia,

    .. math::
        R_g = \sqrt{\frac{2\pi(ABC)^{1/3}N_A}{\text{MW}}}

    Parameters
    ----------
    MW : float
        Molecular weight, [g/mol]
    A : float
        First principle moment of inertia, [kg*m^2]
    B : float
        Second principle moment of inertia, [kg*m^2]
    C : float
        Third principle moment of inertia, [kg*m^2]
    planar : bool
        Whether the molecule is flat or not, [-]

    Returns
    -------
    Rg : float
        Radius of gyration, [m]

    Notes
    -----
    There are many, many quantum chemistry models available which give
    different results.

    Examples
    --------
    Example calcultion in [1]_ for hydrazine (optimized with HF/6-31G model):

    >>> radius_of_gyration(MW=32.00452, planar=False, A=5.692E-47, B=3.367E-46, C=3.681E-46)
    1.50581642e-10

    The same calculation was performed with `psi` and somewhat different parameters obtained

    >>> radius_of_gyration(MW=32.00452, planar=False, A=6.345205205562681e-47, B=3.2663291891213418e-46, C=3.4321304373822523e-46)
    1.507895671e-10

    A planar molecule, bromosilane, has two principle moments of inertia
    in [2]_. They are 2.80700 cm^-1 and 0.14416 cm^-1. These can be converted to
    MHz as follows:

    These can then be converted to units of AMU*Angstrom^2, and from there
    to kg*m^2.

    >>> A, B = 2.80700, 0.14416
    >>> from scipy.constants import atomic_mass, c, angstrom
    >>> A, B = A*c*1e-4, B*c*1e-4 # from cm^-1 to MHz
    >>> A, B = [505379.15/i for i in (A, B)] #  TODO which constants did this conversion factor come from, AMU*Angstrom^2
    >>> A, B = [i*atomic_mass*angstrom**2 for i in (A, B)] # amu*angstrom^2 to kg*m^2
    >>> radius_of_gyration(A=A, B=B, planar=True, MW=111.01, C=0)
    4.8859099776e-11


    Alternatively, doing the conversion all in one:

    >>> A, B = 2.80700, 0.14416
    >>> from scipy.constants import c, h, pi
    >>> A, B = A*c*100, B*c*100 # from cm^-1 to Hz
    >>> A, B = [h/(8*pi**2)/i for i in (A, B)] # from Hz to kg*m^2
    >>> radius_of_gyration(A=A, B=B, planar=True, MW=111.01, C=0)
    4.885909296e-11

    This is also nicely documented on this page: https://cccbdb.nist.gov/convertmomint.asp
    which was unfortunately found by the author after figuring it out the hard way.


    References
    ----------
    .. [1] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       8E. McGraw-Hill Professional, 2007.
    .. [2] Johnson III, Russell D. "NIST 101. Computational Chemistry
       Comparison and Benchmark Database," 1999. https://cccbdb.nist.gov
    '''
    if planar:
        return sqrt(sqrt(A*B)*N_A*1e3/MW)
    else:
        return sqrt(2*pi*(A*B*C)**(1/3)*N_A*1e3/MW)
