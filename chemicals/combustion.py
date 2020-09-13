# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>
Copyright (C) 2020 Yoel Rene Cortes-Pena <yoelcortes@gmail.com>

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


This module contains a series of functions for modeling combustion reactions.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:
    
Combustion Stoichiometry
------------------------
.. autofunction:: chemicals.combustion.combustion_stoichiometry
.. autofunction:: chemicals.combustion.combustion_products_mixture

Heat of Combustion
------------------
.. autofunction:: chemicals.combustion.HHV_stoichiometry
.. autofunction:: chemicals.combustion.HHV_modified_Dulong
.. autofunction:: chemicals.combustion.LHV_from_HHV

Heat of Combustion and Stiochiometry
------------------------------------
.. autofunction:: chemicals.combustion.combustion_data
.. autoclass:: chemicals.combustion.CombustionData

Basic Combustion Spec Solvers
-----------------------------
.. autofunction:: chemicals.combustion.fuel_air_spec_solver
.. autofunction:: chemicals.combustion.combustion_spec_solver
.. autofunction:: chemicals.combustion.air_fuel_ratio_solver

"""

from chemicals.elements import mass_fractions, molecular_weight, simple_formula_parser
from chemicals.utils import property_molar_to_mass, property_mass_to_molar
from fluids.numerics import normalize

__all__ = ('combustion_stoichiometry',
           'CombustionData',
           'combustion_data',
           'HHV_modified_Dulong',
           'HHV_stoichiometry',
           'LHV_from_HHV',
           'combustion_products_mixture',
           'air_fuel_ratio_solver',
           'fuel_air_spec_solver',
           'combustion_spec_solver')

   
def as_atoms(formula):
    if isinstance(formula, str):
        atoms = simple_formula_parser(formula)
    elif isinstance(formula, dict):
        atoms = formula
    else:
        raise ValueError("atoms must be either a string or dictionary, "
                        "not a %s object" %(type(formula).__name__))
    return atoms

DULONG = 'Dulong'
STOICHIOMETRY = 'Stoichiometry'
HHV_methods = (DULONG, STOICHIOMETRY)

combustible_elements = ('C', 'H', 'N', 'O', 'S', 'Br', 'I', 'Cl', 'F', 'P')
combustible_elements_set = frozenset(combustible_elements)

Hf_combustion_chemicals = {
    'H2O': -285825.0,
    'CO2': -393474.0,
    'SO2': -296800.0,
    'Br2': 30880.0,
    'I2': 62417.0,
    'HCl': -92173.0,
    'HF': -272711.0,
    'P4O10': -3009940.0,
    'O2': 0.0,
    'N2': 0.0,
    "Ash": 0.0,
}


incombustible_materials = {
'7440-37-1': 'Ar',
'7782-44-7': 'O2',
'7440-01-9': 'Ne',
'7439-90-9': 'Kr',
'7440-63-3': 'Xe',
'7727-37-9': 'N2',
'124-38-9': 'CO2',
'1314-13-2': 'ZnO2',
'7732-18-5': 'water',
'7789-20-0': 'water(D2)',
'13463-67-7': 'TiO2',
'14762-55-1': 'He3',
'7440-59-7': 'He',
'7782-50-5': 'Cl',
'7446-09-5': 'SO2',
'7726-95-6': 'Br'}

combustion_products_to_CASs = {
'Br2': '7726-95-6',
'CO2': '124-38-9',
'H2O': '7732-18-5',
'HCl': '7647-01-0',
'HF': '7664-39-3',
'I2': '7553-56-2',
'N2': '7727-37-9',
'O2': '7782-44-7',
'P4O10': '16752-60-6',
'SO2': '7446-09-5'}

unreactive_CASs = {
    '124-38-9': 'CO2',
    '16752-60-6': 'P4O10',
    '7446-09-5': 'SO2',
    '7553-56-2': 'I2',
    '7647-01-0': 'HCl',
    '7664-39-3': 'HF',
    '7726-95-6': 'Br2',
    '7727-37-9': 'N2',
    '7732-18-5': 'H2O'}

O2_CAS = '7782-44-7'
H2O_CAS = '7732-18-5'

def combustion_stoichiometry(atoms, MW=None, missing_handling='elemental'):
    r"""Return a dictionary of stoichiometric coefficients of chemical 
    combustion, given a dictionary of a molecule's constituent atoms and their 
    counts. 
    
    This function is based on the combustion of hydrocarbons; the products for 
    some inorganics can be hard to predict, and no special handling is included
    here for them. This reaction is the standard one at standard pressure with
    an excess of oxygen; it does not account for partial combustion or nitrous
    oxides.
    
    Parameters
    ----------
    atoms : dict[str, int]
        Dictionary of atoms and their counts, [-]
    MW : float, optional
        Molecular weight of chemical, used only if `missing_handling` is 'Ash', 
        [g/mol]
    missing_handling : str, optional
        How to handle compounds which do not appear in the stoichiometric
        reaction below. If 'elemental', return those atoms in the monatomic 
        state; if 'ash', converts all missing attoms to 'Ash' in the output at
        a `MW` of 1 g/mol, [-]

    Returns
    -------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of combustion. May inlcude the following 
        keys for complete combustion: 'H2O', 'CO2', 'SO2', 'Br2', 'I2', 'HCl', 
        'HF' 'P4O10'; if `missing_handling` is 'elemental' can include the 
        other elements; if `missing_handling` is 'ash', Ash will be present in 
        the output if the compounds whose reactions are not included here.
        'O2' is always present, with negative values indicating oxygen is
        required. [-]

    Notes
    -----
    The stoichiometry is given by:
        
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O

    .. math::
        k = c + s + \frac{h}{4} + \frac{5P}{4} - \frac{x + f}{4} - \frac{o}{2}

    Also included in the results is the moles of O2 required per mole of
    the mixture of the molecule.
    
    HF and HCl are gaseous products in their standard state. P4O10 is a solid
    in its standard state. Bromine is a liquid as is iodine. Water depends on
    the chosen definition of heating value. The other products are gases.
    
    Atoms not in ['C', 'H', 'N', 'O', 'S', 'Br', 'I', 'Cl', 'F', 'P'] are 
    returned as pure species; i.e. sodium hydroxide produces water and pure
    Na.

    Examples
    --------
    Methane gas burning:
    
    >>> combustion_stoichiometry({'C': 1, 'H':4})
    {'CO2': 1, 'O2': -2.0, 'H2O': 2.0}

    """
    products = {}
    nC, nH, nN, nO, nS, nBr, nI, nCl, nF, nP = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if 'C' in atoms:
        products['CO2'] = nC = atoms['C']
    if 'H' in atoms:
        nH = atoms['H']
    if 'N' in atoms:
        nN = atoms['N']
    if 'O' in atoms:
        nO = atoms['O']
    if 'S' in atoms:
        nS = atoms['S']
    if 'Br' in atoms:
        nBr = atoms['Br']
    if 'I' in atoms:
        nI = atoms['I']
    if 'Cl' in atoms:
        nCl = atoms['Cl']
    if 'F' in atoms:
        nF = atoms['F']
    if 'P' in atoms:
        nP = atoms['P']

    nO2 = -(nC + nS + .25*nH + 1.25*nP - .25*(nCl + nF) - .5*nO)
    nCO2 = nC
    nBr2 = .5*nBr
    nI2 = .5*nI
    nHCl = nCl
    nHF = nF
    nSO2 = nS
    nN2 = .5*nN
    nP4O10 = .25*nP
    nH2O = 0.5*(nH - nCl - nF)
    
    if nO2 != 0.0:
        products['O2'] = nO2
    if nCO2 !=  0.0:
        products['CO2'] = nCO2
    if nBr2 != 0.0:
        products['Br2'] = nBr2
    if nI2 != 0.0:
        products['I2'] = nI2
    if nCl != 0.0:
        products['HCl'] = nHCl
    if nHF != 0.0:
        products['HF'] = nHF
    if nSO2 != 0.0:
        products['SO2'] = nSO2
    if nN2 != 0.0:
        products['N2'] = nN2
    if nP4O10 != 0.0:
        products['P4O10'] = nP4O10
    if nH2O != 0.0:
        products['H2O'] = nH2O
    
    missing_handling = missing_handling.lower() 
    if missing_handling == 'elemental':
        for atom, value in atoms.items():
            if atom not in combustible_elements_set:
                products[atom] = value
    elif missing_handling == 'ash':
        combustion_atoms = {i: atoms.get(i, 0) for i in combustible_elements}        
        MW = MW or molecular_weight(atoms)
        Ash = MW - molecular_weight(combustion_atoms)
        if Ash/MW > 0.0001:
            products['Ash'] = Ash
    else:
        raise ValueError("Allowed values for `missing_handling` are 'elemental' and 'ash'.")
    return products


def combustion_products_mixture(atoms_list, zs, reactivities=None, CASs=None,
                                missing_handling='elemental', 
                                combustion_stoichiometries=None):
    """Calculates the combustion products of a mixture of molecules and their,
    mole fractions; requires a list of dictionaries of each molecule's
    constituent atoms and their counts. Products for non-hydrocarbons may not be
    correct, but are still calculated.

    Parameters
    ----------
    atoms_list : list[dict]
        List of dictionaries of atoms and their counts, [-]
    zs : list[float]
        Mole fractions of each molecule in the mixture, [-]
    reactivities : list[bool]
        Indicators as to whether to combust each molecule, [-]
    CASs : list[str]
        CAS numbers of all compounds; non-reacted products will appear
        in the products indexed by their CAS number, [-]
    missing_handling : str, optional
        How to handle compounds which do not appear in the stoichiometric
        reaction below. If 'elemental', return those atoms in the monatomic 
        state; if 'Ash', converts all missing attoms to 'Ash' in the output at
        a `MW` of 1 g/mol, [-]
    combustion_stoichiometries : list[dict[str, float]]
        List of return values from `combustion_stoichiometry`, can be 
        provided if precomputed [-]

    Returns
    -------
    combustion_producucts : dict
        Dictionary of combustion products and their counts, [-]

    Notes
    -----
    Also included in the results is the moles of O2 required per mole of
    the mixture to be burnt.
    
    Note that if O2 is in the feed, this will be subtracted from the required
    O2 amount.

    HF and HCl are gaseous products in their standard state. P4O10 is a solid
    in its standard state. Bromine is a liquid as is iodine. Water depends on
    the chosen definition of heating value. The other products are gases.
    
    Note that if instead of mole fractions, mole flows are given - the results
    are in terms of mole flows as well!

    Examples
    --------
    Mixture of methane and ethane.

    >>> combustion_products_mixture([{'H': 4, 'C': 1}, {'H': 6, 'C': 2}, {'Ar': 1}, {'C': 15, 'H': 32}],
    ... [.9, .05, .04, .01], reactivities=[True, True, True, False], 
    ... CASs=['74-82-8', '74-84-0', '7440-37-1', '629-62-9'])
    {'CO2': 1.0, 'O2': -1.975, 'H2O': 1.9500000000000002, 'Ar': 0.04, '629-62-9': 0.01}
    """
    # Attempted to use a .copy() on a base dict but that was slower
    products = {}
    has_reactivities = reactivities is not None    
    for i, (atoms, zs_i) in enumerate(zip(atoms_list, zs)):
        if has_reactivities and not reactivities[i]:
            products[CASs[i]] = zs_i
        else:
            ans = (combustion_stoichiometry(atoms, missing_handling=missing_handling)
                   if not combustion_stoichiometries else combustion_stoichiometries[i])
            for key, val in ans.items():
                if key in products:
                    products[key] += val*zs_i
                else:
                    products[key] = val*zs_i
    return products


def combustion_products_to_list(products, CASs):
    zs = [0.0 for i in CASs]
    for product, zi in products.items():
        if product == 'O2_required':
            product = 'O2'
            zi = -zi
            
        if product in CASs:
            zs[CASs.index(product)] = zi
        elif combustion_products_to_CASs[product] in CASs:
            zs[CASs.index(combustion_products_to_CASs[product])] = zi
        else:
            if abs(zi) > 0:
                raise ValueError("Combustion product not in package")
        
    return zs



def is_combustible(CAS, atoms, reactive=True):
    if not reactive:
        return False
    if CAS in unreactive_CASs:
        return False
    elif 'C' in atoms and atoms['C'] > 0:
        return True
    elif 'H' in atoms and atoms['H'] > 0:
        return True
    else:
        return False

def HHV_stoichiometry(stoichiometry, Hf, Hf_chemicals=None):
    r"""
    Return the higher heating value [HHV; in J/mol] based on the 
    theoretical combustion stoichiometry and the heat of formation of
    the chemical.
    
    Parameters
    ----------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of combustion. May inlcude the following 
        keys: 'H2O', 'CO2', 'SO2', 'Br2', 'I2', 'HCl', 'HF' and 'P4O10'.
    Hf : float
        Heat of formation [J/mol].
    Hf_chemicals : dict[str, float]
        Heat of formation of chemicals present in stoichiometry.
    
    Returns
    -------
    HHV : float
        Higher heating value [J/mol].
    
    Notes
    -----
    The combustion reaction is based on the following equation:
    
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O
        
    .. math::
        k = c + s + \frac{h}{4} + \frac{5P}{4} - \frac{x + f}{4} - \frac{o}{2}
    
    The HHV is calculated as the heat of reaction.
    
    Examples
    --------
    Burning methane gas:
        
    >>> HHV_stoichiometry({'O2': -2.0, 'CO2': 1, 'H2O': 2.0}, -74520.0)
    -890604.0
    
    """
    Hfs = Hf_chemicals or Hf_combustion_chemicals
    return sum([Hfs[i] * j for i, j in stoichiometry.items()]) - Hf

def HHV_modified_Dulong(mass_fractions):
    r"""
    Return higher heating value [HHV; in J/g] based on the modified 
    Dulong's equation [1]_.
    
    Parameters
    ----------
    mass_fractions : dict[str, float]
        Dictionary of atomic mass fractions [-].
    
    Returns
    -------
    HHV : float
        Higher heating value [J/mol].
    
    Notes
    -----
    The heat of combustion in J/mol is given by Dulong's equation [1]_:
    
    .. math:: 
        Hc (J/mol) = MW \cdot (338C + 1428(H - O/8)+ 95S)
    
    This equation is only good for <10 wt. % Oxygen content. Variables C, H, O,
    and S are atom weight fractions.
    
    Examples
    --------
    Dry bituminous coal:
    
    >>> HHV_modified_Dulong({'C': 0.716, 'H': 0.054, 'S': 0.016, 'N': 0.016, 'O': 0.093, 'Ash': 0.105})
    -304.0395
    
    References
    ----------
    .. [1] Green, D. W. Waste management. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018
    
    """
    C = mass_fractions.get('C', 0.)
    H = mass_fractions.get('H', 0.)
    O = mass_fractions.get('O', 0.)
    S = mass_fractions.get('S', 0.)
    if O > 0.105:
        raise ValueError("Dulong's formula is only valid at 10 wt. %% Oxygen "
                         "or less (%s given)" %(O))
    return - (338.*C  + 1428.*(H - O/8.)+ 95.*S)    

def LHV_from_HHV(HHV, N_H2O):
    r"""
    Return the lower heating value [LHV; in J/mol] of a chemical given
    the higher heating value [HHV; in J/mol] and the number of water
    molecules formed per molecule burned.
    
    Parameters
    ----------
    HHV : float
        Higher heating value [J/mol].
    N_H2O : int
        Number of water molecules produced [-].
    
    Returns
    -------
    LHV : float
        Lower heating value [J/mol].
    
    Notes
    -----
    The LHV is calculated as follows:
        
    .. math::
        LHV = HHV + H_{vap} \cdot H_2O
        
    .. math::
        H_{vap} = 44011.496 \frac{J}{mol H_2O}
        
    .. math::
        H_2O = \frac{mol H_2O}{mol}
    
    Examples
    --------
    Methanol lower heat of combustion:
    
    >>> LHV_from_HHV(-726024.0, 2)
    -638001.008
    
    """
    return HHV + 44011.496 * N_H2O

def combustion_data(formula=None, stoichiometry=None, Hf=None, MW=None,
                    method=None, missing_handling='ash'):
    r"""
    Return a CombustionData object (a named tuple) that contains the stoichiometry 
    coefficients of the reactants and products, the lower and higher 
    heating values [LHV, HHV; in J/mol], the heat of formation [Hf; in J/mol],
    and the molecular weight [MW; in g/mol].

    Parameters
    ----------
    formula : str, or dict[str, float], optional
        Chemical formula as a string or a dictionary of atoms and their counts.
    stoichiometry : dict[str, float], optional
        Stoichiometry of combustion reaction.
    Hf : float, optional
        Heat of formation of given chemical [J/mol].
        Required if method is "Stoichiometry".
    MW : float, optional
        Molecular weight of chemical [g/mol].
    method : "Stoichiometry" or "Dulong", optional
        Method to estimate LHV and HHV.    
    missing_handling : str, optional
        How to handle compounds which do not appear in the stoichiometric
        reaction below. If 'elemental', return those atoms in the monatomic 
        state; if 'Ash', converts all missing attoms to 'Ash' in the output at
        a `MW` of 1 g/mol, [-]
        
    Returns
    -------
    combustion_data : :class:`~chemicals.combustion.CombustionData`
        A combustion data object with the stoichiometric coefficients of 
        combustion, higher heating value, heat of formation, and molecular 
        weight as attributes named stoichiomery, HHV, Hf, and MW, respectively.
    
    Notes
    -----
    The combustion reaction is based on the following equation:
    
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O
        
    .. math::
        k = c + s + \frac{h}{4} + \frac{5P}{4} - \frac{x + f}{4} - \frac{o}{2}

    If the method is "Stoichiometry", the HHV is found using 
    through an energy balance on the reaction (i.e. heat of reaction).
    If the method is "Dulong", Dulong's equation is used [1]_:

    .. math:: 
        Hc (J/mol) = MW \cdot (338C + 1428(H - O/8)+ 95S)
        
    The LHV is calculated as follows:
        
    .. math::
        LHV = HHV + H_{vap} \cdot H_2O
        
    .. math::
        H_{vap} = 44011.496 \frac{J}{mol H_2O}
        
    .. math::
        H_2O = \frac{mol H_2O}{mol}
        
    Examples
    --------
    Liquid methanol burning:

    >>> combustion_data({'H': 4, 'C': 1, 'O': 1}, Hf=-239100)
    CombustionData(stoichiometry={'CO2': 1, 'O2': -1.5, 'H2O': 2.0}, HHV=-726024.0, Hf=-239100, MW=32.04186)
    
    References
    ----------
    .. [1] Green, D. W. Waste management. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018
    
    """
    if formula:
        if stoichiometry:
            raise ValueError("must specify either `formula` or `stoichiometry`, not both")
        atoms = as_atoms(formula)
    if not stoichiometry:
        try: 
            stoichiometry = combustion_stoichiometry(atoms, MW, missing_handling)
        except NameError:
            raise ValueError("must specify either `formula` or `stoichiometry`, none specified")
    if MW is None:
        MW = molecular_weight(atoms)
    if method:
        method = method.capitalize()
    else:
        method = 'Dulong' if Hf is None else 'Stoichiometry'
    if method == DULONG:
        HHV = MW * HHV_modified_Dulong(mass_fractions(atoms))
        if Hf: raise ValueError("cannot specify Hf if method is 'Dulong'")
        Hf = HHV - HHV_stoichiometry(stoichiometry, 0)
    elif method == STOICHIOMETRY:
        if Hf is None: raise ValueError("must specify Hf if method is 'Stoichiometry'")
        HHV = HHV_stoichiometry(stoichiometry, Hf)
    else:
        raise ValueError("method must be either 'Stoichiometric' or 'Dulong', "
                         "not %s" %(method))
    return CombustionData(stoichiometry, HHV, Hf, MW)


class CombustionData(object):
    r"""
    Return a CombustionData object (a named tuple) that contains the stoichiometry 
    coefficients of the reactants and products, the lower and higher 
    heating values [LHV, HHV; in J/mol], the heat of formation [Hf; in J/mol],
    and the molecular weight [MW; in g/mol].
    
    Parameters
    ----------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of the reactants and products.
    HHV : float
        Higher heating value [J/mol].
    Hf : float
        Heat of formation [J/mol].
    MW : float
        Molecular weight [g/mol].
    
    """
    
    def __init__(self, stoichiometry, HHV, Hf, MW):
        self.stoichiometry = stoichiometry
        self.HHV = HHV
        self.Hf = Hf
        self.MW = MW

    @property
    def LHV(self):
        """Lower heating value [LHV; in J/mol]"""
        return LHV_from_HHV(self.HHV, self.stoichiometry.get('H2O', 0.))

    def __repr__(self):
        return 'CombustionData(stoichiometry=%s, HHV=%s, Hf=%s, MW=%s)' % (
                    self.stoichiometry, self.HHV, self.Hf, self.MW)


def air_fuel_ratio_solver(ratio, Vm_air, Vm_fuel, MW_air, MW_fuel,
                          n_air=None, n_fuel=None,
                          basis='mass'):
    """Calculates molar flow rate of air or fuel from the other, using a
    specified air-fuel ratio. Supports 'mole', 'mass', and 'volume'.

    bases for the ratio variable. The ratio must be of the same units - 
    i.e. kg/kg instead of lb/kg.
    
    The mole, mass, and volume air-fuel ratios are calculated in the process
    and returned as well.

    Parameters
    ----------
    ratio : float
        Air-fuel ratio, in the specified `basis`, [-]
    Vm_air : float
        Molar volume of air, [m^3/mol]
    Vm_fuel : float
        Molar volume of fuel, [m^3/mol]
    MW_air : float
        Molecular weight of air, [g/mol]
    MW_fuel : float
        Molecular weight of fuel, [g/mol]
    n_air : float, optional
        Molar flow rate of air, [mol/s]
    n_fuel : float, optional
        Molar flow rate of fuel, [mol/s]
    basis : str, optional
        One of 'mass', 'mole', or 'volume', [-]
        
    Returns
    -------
    n_air : float
        Molar flow rate of air, [mol/s]
    n_fuel : float
        Molar flow rate of fuel, [mol/s]
    mole_ratio : float
        Air-fuel mole ratio, [-]
    mass_ratio : float
        Air-fuel mass ratio, [-]
    volume_ratio : float
        Air-fuel volume ratio, [-]

    Notes
    -----
    The function works so long as the flow rates, molar volumes, and molecular
    weights are in a consistent basis.
    
    The function may also be used to obtain the other ratios, even if both 
    flow rates are known.
    
    Be careful to use standard volumes if the ratio known is at standard
    conditions! 
    
    This function has no provision for mixed units like mass/mole or 
    volume/mass.
    
    Examples
    --------
    >>> Vm_air = 0.024936627188566596
    >>> Vm_fuel = 0.024880983160354486
    >>> MW_air = 28.850334
    >>> MW_fuel = 17.86651
    >>> n_fuel = 5.0
    >>> n_air = 25.0
    >>> air_fuel_ratio_solver(ratio=5.0, Vm_air=Vm_air, Vm_fuel=Vm_fuel,
    ... MW_air=MW_air, MW_fuel=MW_fuel, n_air=n_air,
    ... n_fuel=n_fuel, basis='mole')
    (25.0, 5.0, 5.0, 8.073858296891782, 5.011182039683378)
    """
    if basis == 'mole':
        if n_air is not None and n_fuel is None:
            n_fuel = n_air/ratio
        elif n_fuel is not None and n_air is None:
            n_air = n_fuel*ratio
    elif basis == 'mass':
        if n_air is not None and n_fuel is None:
            m_air = property_mass_to_molar(n_air, MW_air)
            m_fuel = m_air/ratio
            n_fuel = property_molar_to_mass(m_fuel, MW_fuel)
        elif n_fuel is not None and n_air is None:
            m_fuel = property_mass_to_molar(n_fuel, MW_fuel)
            m_air = m_fuel*ratio
            n_air = property_molar_to_mass(m_air, MW_air)
    elif basis == 'volume':
        if n_air is not None and n_fuel is None:
            V_air = n_air*Vm_air
            V_fuel = V_air/ratio
            n_fuel = V_fuel/Vm_fuel
        elif n_fuel is not None and n_air is None:
            V_fuel = n_fuel*Vm_fuel
            V_air = V_fuel*ratio
            n_air = V_air/Vm_air
    if n_air is None or n_fuel is None:
        raise ValueError("Could not convert")
    mole_ratio = n_air/n_fuel
    mass_ratio, volume_ratio = MW_air/MW_fuel*mole_ratio, Vm_air/Vm_fuel*mole_ratio
    return n_air, n_fuel, mole_ratio, mass_ratio, volume_ratio


def fuel_air_spec_solver(zs_air, zs_fuel, CASs, atomss, n_fuel=None, 
                         n_air=None, n_out=None,
                         O2_excess=None, frac_out_O2=None,
                         frac_out_O2_dry=None, ratio=None,
                         Vm_air=None, Vm_fuel=None, MW_air=None, MW_fuel=None,
                         ratio_basis='mass', reactivities=None,
                         combustion_stoichiometries=None):
    """Solves the system of equations describing a flow of air mixing with a
    flow of combustibles and burning completely. All calculated variables are
    returned as a dictionary.

    Supports solving with any 2 of the extensive variables, or one extensive
    and one intensive variable:
    
    Extensive variables:
        
    * `n_air`
    * `n_fuel`
    * `n_out`
    
    Intensive variables:
        
    * `O2_excess`
    * `frac_out_O2`
    * `frac_out_O2_dry`
    * `ratio`
    
    The variables `Vm_air`, `Vm_fuel`, `MW_air`, and `MW_fuel` are only
    required when an air-fuel ratio is given. Howver, the ratios cannot be
    calculated for the other solve options without them.
    
    Parameters
    ----------
    zs_air : list[float]
        Mole fractions of the air; most not contain any combustibles, [-]
    zs_fuel : list[float]
        Mole fractions of the fuel; can contain inerts and/or oxygen as well,
        [-]
    CASs : list[str]
        CAS numbers of all compounds, [-]
    atomss : list[dict[float]]
        List of dictionaries of elements and their counts for all molecules in
        the mixtures, [-]
    n_fuel : float, optional
        Flow rate of fuel, [mol/s]
    n_air : float, optional
        Flow rate of air, [mol/s]
    n_out : float, optional
        Flow rate of combustion products, remaining oxygen, and inerts, [mol/s]
    O2_excess : float, optional
        The excess oxygen coming out; (O2 in)/(O2 required) - 1, [-]
    frac_out_O2 : float, optional
        The mole fraction of oxygen out, [-]
    frac_out_O2_dry : float, optional
        The mole fraction of oxygen out on a dry basis, [-]
    ratio : float, optional
        Air-fuel ratio, in the specified `basis`, [-]
    Vm_air : float, optional
        Molar volume of air, [m^3/mol]
    Vm_fuel : float, optional
        Molar volume of fuel, [m^3/mol]
    MW_air : float, optional
        Molecular weight of air, [g/mol]
    MW_fuel : float, optional
        Molecular weight of fuel, [g/mol]
    ratio_basis : str, optional
        One of 'mass', 'mole', or 'volume', [-]
    reactivities : list[bool], optional
        Optional list which can be used to mark otherwise combustible 
        compounds as incombustible and which will leave unreacted, [-]
    combustion_stoichiometries : list[dict[str, float]]
        List of return values from `combustion_stoichiometry`, can be 
        provided if precomputed [-]
        
    Returns
    -------
    results : dict
        * n_fuel : Flow rate of fuel, [mol/s]
        * n_air : Flow rate of air, [mol/s]
        * n_out : Flow rate of combustion products, remaining oxygen, and 
          inerts, [mol/s]
        * O2_excess : The excess oxygen coming out; (O2 in)/(O2 required) - 1,
          [-]
        * frac_out_O2 : The mole fraction of oxygen out, [-]
        * frac_out_O2_dry : The mole fraction of oxygen out on a dry basis, [-]
        * mole_ratio : Air-fuel mole ratio, [-]
        * mass_ratio : Air-fuel mass ratio, [-]
        * volume_ratio : Air-fuel volume ratio, [-]
        * ns_out : Mole flow rates out, [mol/s] 
        * zs_out : Mole fractions out, [-]

    Notes
    -----
    Combustion products themselves cannot be set as unreactive.
    
    The function works so long as the flow rates, molar volumes, and molecular
    weights are in a consistent basis.
    
    The function may also be used to obtain the other ratios, even if both 
    flow rates are known.
    
    Be careful to use standard volumes if the ratio known is at standard
    conditions! 
    
    Examples
    --------
    >>> zs_air = [0.79, 0.205, 0, 0, 0, 0.0045, 0.0005]
    >>> zs_fuel = [0.025, 0.025, 0.85, 0.07, 0.029, 0.0005, 0.0005]
    >>> CASs = ['7727-37-9', '7782-44-7', '74-82-8', '74-84-0', '74-98-6', '7732-18-5', '124-38-9']
    >>> atomss = [{'N': 2}, {'O': 2}, {'H': 4, 'C': 1}, {'H': 6, 'C': 2}, {'H': 8, 'C': 3}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}]
    >>> ans = fuel_air_spec_solver(zs_air=zs_air, zs_fuel=zs_fuel, CASs=CASs, atomss=atomss, n_fuel=5.0, O2_excess=0.3, Vm_air=0.02493, Vm_fuel=0.02488, MW_air=28.79341351, MW_fuel=18.55158039)
    >>> [round(i, 5) for i in ans['ns_out']]
    [51.99524, 3.135, 0.0, 0.0, 0.0, 10.42796, 5.42033]
    >>> [round(i, 5) for i in ans['zs_out']]
    [0.73255, 0.04417, 0.0, 0.0, 0.0, 0.14692, 0.07637]
    >>> ans['frac_out_O2'], ans['frac_out_O2_dry']
    (0.04416828172034148, 0.051774902132807)
    >>> ans['mole_ratio'], ans['mass_ratio'], ans['volume_ratio']
    (13.131707317073175, 20.381372957130615, 13.15809740412517)
    >>> ans['n_air']
    65.65853658536588
    """
    TRACE_FRACTION_IN_AIR = 1e-10
    # Only one path to get n_air, n_fuel should ever be followed. Calculate all the
    # extra information redundantly at the end!

    # Handle combustibles in the air by burning them right away,
    # and working with that n_air.
    
    N = len(CASs)
    cmps = range(N)
    
    if reactivities is None:
        reactivities = [True for i in zs_air]
    combustibilities = [is_combustible(CASs[i], atomss[i], reactivities[i]) for i in cmps]
    
    for i in combustibilities:
        if zs_air[i] > TRACE_FRACTION_IN_AIR:
            pass
        
        
    O2_index = CASs.index(O2_CAS)
    H2O_index = CASs.index(H2O_CAS)
    
    z_air_O2 = zs_air[O2_index]
    z_air_H20 = zs_air[H2O_index]
    z_fuel_O2 = zs_fuel[O2_index]
    
    if ratio is not None and (n_air is None or n_fuel is None):
        n_air, n_fuel = air_fuel_ratio_solver(ratio, Vm_air, Vm_fuel, MW_air, MW_fuel, n_air=n_air,
                                              n_fuel=n_fuel, basis=ratio_basis)[0:2]
    
    # Given O2 excess and either air or fuel flow rate, can solve directly for the other
    if O2_excess is not None and (n_fuel is None or n_air is None):
        if n_fuel is not None:
            comb_ans = combustion_products_mixture(atomss, zs_fuel, reactivities=reactivities, CASs=CASs,
                                                   combustion_stoichiometries=combustion_stoichiometries)
            n_O2_required = -comb_ans['O2']*n_fuel + z_fuel_O2*n_fuel
            N_O2_required_air = n_O2_required*(1.0 + O2_excess)- z_fuel_O2*n_fuel
            n_air = N_O2_required_air/zs_air[O2_index]
        elif n_air is not None:
            '''from sympy import *
            O2, n_fuel, O2_per_mole_fuel, n_O2_air, z_fuel_O2 = symbols('O2, n_fuel, O2_per_mole_fuel, n_O2_air, z_fuel_O2')
            Eq1 = Eq(O2/n_fuel, O2_per_mole_fuel)
            Eq2 = Eq(O2, n_O2_air + z_fuel_O2*n_fuel)
            solve([Eq1, Eq2], [n_fuel, O2])'''
            comb_ans = combustion_products_mixture(atomss, zs_fuel, reactivities=reactivities, CASs=CASs,
                                                   combustion_stoichiometries=combustion_stoichiometries)
            O2_per_mole_fuel = (-comb_ans['O2'] + zs_fuel[O2_index])*(1.0 + O2_excess)
            n_O2_air = zs_air[O2_index]*n_air
            n_fuel = n_O2_air/(O2_per_mole_fuel - z_fuel_O2)
        elif n_out is not None:
            '''from sympy import *
            n_air, n_fuel, n_out, n_delta, O2_coeff, O2_excess, z_air_O2, z_fuel_O2 = symbols('n_air, n_fuel, n_out, n_delta, O2_coeff, O2_excess, z_air_O2, z_fuel_O2')
            n_O2_in = n_air*z_air_O2 + z_fuel_O2*n_fuel
            Eq1 = Eq(O2_excess, n_O2_in/(n_fuel*O2_coeff) - 1)
            Eq2 = Eq(n_out, n_air + (n_delta)*n_fuel + n_fuel)
            solve([Eq1, Eq2], [n_fuel, n_air])'''
            comb_ans = combustion_products_mixture(atomss, zs_fuel, reactivities=reactivities, CASs=CASs,
                                                   combustion_stoichiometries=combustion_stoichiometries)
            stoic_comb_products = combustion_products_to_list(comb_ans, CASs)
            O2_coeff = O2_burnt_n_fuel = -comb_ans['O2'] + z_fuel_O2
            n_delta = 1
            for z_new, z_old, CAS in zip(stoic_comb_products, zs_fuel, CASs):
                if CAS != O2_CAS:
                    n_delta += z_new - z_old
            n_delta -= -comb_ans['O2'] + z_fuel_O2

            n_fuel = n_out*z_air_O2/(O2_coeff*O2_excess + O2_coeff + n_delta*z_air_O2 - z_fuel_O2)
            n_air = n_out*(O2_coeff*O2_excess + O2_coeff - z_fuel_O2)/(O2_coeff*O2_excess + O2_coeff + n_delta*z_air_O2 - z_fuel_O2)
            
#            n_fuel = n_out*z_air_O2/(O2_burnt_n_fuel*O2_excess + O2_burnt_n_fuel + n_delta*z_air_O2)
#            n_air = O2_burnt_n_fuel*n_out*(O2_excess + 1)/(O2_burnt_n_fuel*O2_excess + O2_burnt_n_fuel + n_delta*z_air_O2)          
            
            
            
    # Solvers for frac_out_O2 and dry basis
    if (frac_out_O2 is not None or frac_out_O2_dry is not None) and (n_fuel is None or n_air is None):
        if n_fuel is not None:
            ns_fuel = [zi*n_fuel for zi in zs_fuel]
            comb_ans = combustion_products_mixture(atomss, ns_fuel, reactivities=reactivities, CASs=CASs,
                                                   combustion_stoichiometries=combustion_stoichiometries)
            n_O2_stoic = -comb_ans['O2']
            if n_O2_stoic < 0.0:
                raise ValueError("Cannot meet air spec - insufficient air for full combustion")
            # when burning stoichiometrically, how many moles out?
            stoic_comb_products = combustion_products_to_list(comb_ans, CASs)
            stoic_comb_products[O2_index] = 0 # Set the O2 to zero as it is negative
            n_fixed_products = sum(stoic_comb_products)

            # Following two equations solved using SymPy
            if frac_out_O2 is not None:
                '''from sympy import *
                from sympy.abc import *
                frac_goal, n_air, n_fixed, z_O2, O2_burnt = symbols('frac_goal, n_air, n_fixed, z_O2, O2_burnt')
                solve(Eq(frac_goal, (z_O2*n_air-O2_burnt)/((n_fixed+n_air)-O2_burnt)), n_air)
                '''
                n_air = (n_O2_stoic*frac_out_O2 - n_O2_stoic - frac_out_O2*n_fixed_products)/(frac_out_O2 - z_air_O2)
            elif frac_out_O2_dry is not None:
                '''from sympy import *
                from sympy.abc import *
                frac_goal, n_air, n_fixed, z_O2, n_O2_stoic, z_H20 = symbols('frac_goal, n_air, n_fixed, z_O2, n_O2_stoic, z_H20')
                solve(Eq(frac_goal, (z_O2*n_air-n_O2_stoic)/((n_fixed+n_air)-n_O2_stoic -n_air*z_H20)), n_air)
                '''
                n_fixed_products -= stoic_comb_products[H2O_index]
                n_air = ((-frac_out_O2_dry*n_O2_stoic + frac_out_O2_dry*n_fixed_products + n_O2_stoic)
                         /(frac_out_O2_dry*z_air_H20 - frac_out_O2_dry + z_air_O2))
        elif n_air is not None or n_out is not None:
            comb_ans = combustion_products_mixture(atomss, zs_fuel, reactivities=reactivities, CASs=CASs,
                                                   combustion_stoichiometries=combustion_stoichiometries)
            O2_burnt_n_fuel = -comb_ans['O2']
            H2O_from_fuel_n = comb_ans['H2O']
            n_delta = 0.0
            for k, v in comb_ans.items():
                n_delta += v
            
            if frac_out_O2 is not None and n_air is not None:
                '''from sympy import *
                from sympy.abc import *
                frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff = symbols('frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff')
                n_unburnt_O2 = n_air*z_O2 - n_fuel*O2_coeff
                n_out = n_air + coeff*n_fuel
                solve(Eq(frac_goal, n_unburnt_O2/n_out), n_fuel)
                '''
                n_fuel = n_air*(-frac_out_O2 + z_air_O2)/(O2_burnt_n_fuel + n_delta*frac_out_O2)
            elif frac_out_O2_dry is not None and n_air is not None:
                '''from sympy import *
                from sympy.abc import *
                frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, H2O_coeff, z_H2O = symbols('frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, H2O_coeff, z_H2O')
                n_unburnt_O2 = n_air*z_O2 - n_fuel*O2_coeff
                n_out = n_air + coeff*n_fuel
                n_H2O = n_air*z_H2O + H2O_coeff*n_fuel
                solve(Eq(frac_goal, n_unburnt_O2/(n_out-n_H2O)), n_fuel)'''
                n_fuel = n_air*(frac_out_O2_dry*z_air_H20 - frac_out_O2_dry + z_air_O2)/(-H2O_from_fuel_n*frac_out_O2_dry + O2_burnt_n_fuel + n_delta*frac_out_O2_dry)
            elif frac_out_O2 is not None and n_out is not None:
                '''from sympy import *
                frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, n_out = symbols('frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, n_out')
                n_unburnt_O2 = n_air*z_O2 - n_fuel*O2_coeff
                Eq1 = Eq(frac_goal, n_unburnt_O2/n_out)
                Eq2 = Eq(n_out, n_air + coeff*n_fuel)
                solve([Eq1, Eq2], [n_fuel, n_air])'''
                n_fuel = -n_out*(frac_out_O2 - z_air_O2)/(O2_burnt_n_fuel + n_delta*z_air_O2)
                n_air = n_out*(O2_burnt_n_fuel + n_delta*frac_out_O2)/(O2_burnt_n_fuel + n_delta*z_air_O2)
            elif frac_out_O2_dry is not None and n_out is not None:
                '''from sympy import *
                frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, n_out, z_H2O, H2O_coeff = symbols('frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, n_out, z_H2O, H2O_coeff')
                n_unburnt_O2 = n_air*z_O2 - n_fuel*O2_coeff
                n_H2O = n_air*z_H2O + H2O_coeff*n_fuel
                Eq1 = Eq(frac_goal, n_unburnt_O2/(n_out-n_H2O))
                Eq2 = Eq(n_out, n_air + coeff*n_fuel)
                solve([Eq1, Eq2], [n_fuel, n_air])'''
                n_fuel = n_out*(frac_out_O2_dry*z_air_H20 - frac_out_O2_dry + z_air_O2)/(-H2O_from_fuel_n*frac_out_O2_dry + O2_burnt_n_fuel + n_delta*frac_out_O2_dry*z_air_H20 + n_delta*z_air_O2)
                n_air = n_out*(-H2O_from_fuel_n*frac_out_O2_dry + O2_burnt_n_fuel + n_delta*frac_out_O2_dry)/(-H2O_from_fuel_n*frac_out_O2_dry + O2_burnt_n_fuel + n_delta*(frac_out_O2_dry*z_air_H20 + z_air_O2))

    # Case of two fuels known - one n out
    if n_out is not None and (n_fuel is None or n_air is None):
        '''from sympy import *
        from sympy.abc import *
        frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, n_out = symbols('frac_goal, n_air, z_O2, n_fuel, coeff, O2_coeff, n_out')
        solve(Eq(n_out, n_air + coeff*n_fuel), n_fuel)
        solve(Eq(n_out, n_air + coeff*n_fuel), n_air)'''
        comb_ans = combustion_products_mixture(atomss, zs_fuel, reactivities=reactivities, CASs=CASs,
                                               combustion_stoichiometries=combustion_stoichiometries)
        n_delta = 0.0
        for k, v in comb_ans.items():
            n_delta += v
                
        if n_fuel is not None:
            n_air = -n_delta*n_fuel + n_out
        elif n_air is not None:
            n_fuel = (-n_air + n_out)/n_delta
            
        
    
    
    # Compute all other properties from the air and fuel flow rate
    results = {'n_fuel': n_fuel, 'n_air': n_air}
    if n_fuel is not None and n_air is not None:
        ns_to_combust = []
        for zi_air, zi_fuel in zip(zs_air, zs_fuel):
            ns_to_combust.append(zi_air*n_air + zi_fuel*n_fuel)
        
        comb_ans = combustion_products_mixture(atomss, ns_to_combust, reactivities=reactivities, CASs=CASs,
                                               combustion_stoichiometries=combustion_stoichiometries)
        ns_out = combustion_products_to_list(comb_ans, CASs)

#        comb_fuel_only = combustion_products_mixture(atomss, [n_fuel*zi for zi in zs_fuel], reactivities=reactivities, CASs=CASs,
#                                                     combustion_stoichiometries=combustion_stoichiometries)
        
        n_out = sum(ns_out)
        zs_out = normalize(ns_out)
        frac_out_O2 = zs_out[O2_index]
        frac_out_O2_dry = frac_out_O2/(1.0 - zs_out[H2O_index])
        
        results['ns_out'] = ns_out
        results['zs_out'] = zs_out
        results['n_out'] = n_out
        results['frac_out_O2'] = frac_out_O2
        results['frac_out_O2_dry'] = frac_out_O2_dry
        
        
        O2_in = n_air*zs_air[O2_index] + n_fuel*zs_fuel[O2_index]
        O2_demand = O2_in - results['ns_out'][O2_index]
        results['O2_excess'] = O2_in/O2_demand - 1.0
#        results['O2_excess'] = n_air*z_air_O2/(comb_fuel_only['O2_required']) - 1
        
        ratios = (None, None, None)
        if Vm_air is not None and Vm_fuel is not None and MW_air is not None and MW_fuel is not None:
            ratios = air_fuel_ratio_solver(None, Vm_air, Vm_fuel, MW_air, MW_fuel, 
                                           n_air=n_air, n_fuel=n_fuel)[2:]
        results['mole_ratio'], results['mass_ratio'], results['volume_ratio'] = ratios
        
    return results


def fuel_air_third_spec_solver(zs_air, zs_fuel, zs_third, CASs, atomss, n_third,
                           n_fuel=None, n_air=None, n_out=None,
                           O2_excess=None, frac_out_O2=None,
                           frac_out_O2_dry=None, ratio=None,
                           Vm_air=None, Vm_fuel=None, Vm_third=None, 
                           MW_air=None, MW_fuel=None, MW_third=None,
                           ratio_basis='mass', reactivities=None,
                           combustion_stoichiometries=None):
    # I believe will always require 1 intensive spec (or flow of )
    # To begin - exclude n_out spec? Should be possible to solve for it/with it though.
    
    common_specs = {'zs_air': zs_air, 'CASs': CASs, 'atomss': atomss,
                    'n_air': n_air, 'O2_excess': O2_excess, 'frac_out_O2': frac_out_O2,
                    'frac_out_O2_dry': frac_out_O2_dry, 'Vm_air': Vm_air, 'MW_air': MW_air,
                    'ratio': ratio, 'ratio_basis': ratio_basis,
                    'reactivities': reactivities}
    
    O2_index = CASs.index(O2_CAS)
    H2O_index = CASs.index(H2O_CAS)
    z_air_O2 = zs_air[O2_index]
    z_fuel_O2 = zs_fuel[O2_index]

    Vm_mix, MW_mix = None, None

    def fix_ratios(results):
        if n_fuel is not None:
            if Vm_fuel is not None and Vm_third is not None:
                Vm_mix = (Vm_fuel*n_fuel + Vm_third*n_third)/(n_fuel+n_third)
            if MW_fuel is not None and MW_third is not None:
                MW_mix = (MW_fuel*n_fuel + MW_third*n_third)/(n_fuel+n_third)
        # Needs n_fuel and n_air to be in outer namespace
        ratios = (None, None, None)
        if Vm_air is not None and Vm_mix is not None and MW_air is not None and MW_mix is not None:
            ratios = air_fuel_ratio_solver(None, Vm_air, Vm_mix, MW_air, MW_mix, 
                                           n_air=n_air, n_fuel=n_fuel+n_third)[2:]
        results['mole_ratio'], results['mass_ratio'], results['volume_ratio'] = ratios
            
    if n_fuel is not None :
        if Vm_fuel is not None and Vm_third is not None:
            Vm_mix = (Vm_fuel*n_fuel + Vm_third*n_third)/(n_fuel + n_third)
        if MW_fuel is not None and MW_third is not None:
            MW_mix = (MW_fuel*n_fuel + MW_third*n_third)/(n_fuel + n_third)

        # Combine the two fuels for one burn
        n_fuel_mix = n_fuel + n_third
        ns_fuel_mix = [zi*n_fuel + zj*n_third for zi, zj in zip(zs_fuel, zs_third)]
        zs_fuel_mix = normalize(ns_fuel_mix)
        
        mix_burn = fuel_air_spec_solver(zs_fuel=zs_fuel_mix, n_fuel=n_fuel_mix, 
                                        Vm_fuel=Vm_mix, MW_fuel=MW_mix,
                                        n_out=n_out, combustion_stoichiometries=combustion_stoichiometries,
                                       **common_specs)
        mix_burn['n_fuel'] -= n_third
        return mix_burn
    if n_air is not None and n_fuel is None:
        O2_in_third = n_third*zs_third[O2_index]
        O2_in_orig = n_air*z_air_O2 + O2_in_third
        
        third_burn = fuel_air_spec_solver(zs_fuel=zs_third, n_fuel=n_third,
                                          Vm_fuel=Vm_third, MW_fuel=MW_third,
                                          combustion_stoichiometries=combustion_stoichiometries,
                                          **common_specs)
        n_air2 = third_burn['n_out']
        zs_air2 = third_burn['zs_out']

        O2_demand_third = O2_in_orig - third_burn['ns_out'][O2_index]
        
        extra_specs = common_specs.copy()
        extra_specs['n_air'] = n_air2
        extra_specs['zs_air'] = zs_air2
        if O2_excess is not None:
            '''from sympy import *
            O2_excess, excess_sub, O2_air, O2_in_third, O2_fixed, z_fuel_O2, n_fuel, O2_coeff, O2_demand_third = symbols('O2_excess, excess_sub, O2_air, O2_in_third, O2_fixed, z_fuel_O2, n_fuel, O2_coeff, O2_demand_third')
            Eq1 = Eq(O2_excess, (O2_air + O2_in_third + z_fuel_O2*n_fuel)/(O2_demand_third + n_fuel*O2_coeff) - 1)
            Eq2 = Eq(excess_sub, (O2_fixed +z_fuel_O2*n_fuel)/(n_fuel*O2_coeff) - 1)
            solve([Eq1, Eq2], [excess_sub, n_fuel])'''
            comb_ans = combustion_products_mixture(atomss, zs_fuel, reactivities=reactivities, CASs=CASs,
                                                   combustion_stoichiometries=combustion_stoichiometries)
            O2_burnt_n_fuel = O2_coeff = -comb_ans['O2'] + z_fuel_O2
            O2_air = z_air_O2*n_air
            # Simpler expression for n_fuel than fake O2 excess
            n_fuel = (O2_air - O2_demand_third*O2_excess - O2_demand_third + O2_in_third)/(O2_coeff*O2_excess + O2_coeff - z_fuel_O2)

        elif ratio is not None:
            if ratio_basis == 'mole':
                n_fuel = n_air/ratio - n_third
            elif ratio_basis == 'mass':
                m_air, m_third = n_air*MW_air*1e-3, n_third*MW_third*1e-3
                m_fuel = m_air/ratio - m_third
                n_fuel = m_fuel/(MW_fuel*1e-3)
            elif ratio_basis == 'volume':
                Q_air, Q_third = n_air*Vm_air, n_third*Vm_third
                Q_fuel = Q_air/ratio - Q_third
                n_fuel = Q_fuel/Vm_fuel
        
        
        fuel_burn = fuel_air_spec_solver(zs_fuel=zs_fuel, n_fuel=n_fuel,
                                         Vm_fuel=Vm_fuel, MW_fuel=MW_fuel,
                                         n_out=n_out, combustion_stoichiometries=combustion_stoichiometries,
                                         **extra_specs)
        
        O2_in_fuel = fuel_burn['n_fuel']*zs_fuel[O2_index]
        O2_demand_fuel = (n_air2*zs_air2[O2_index] + O2_in_fuel)/(fuel_burn['O2_excess'] + 1)
        

        
        # denominator and numerator are correct now, do not change vars used..
        fuel_burn['O2_excess'] = (O2_in_orig + O2_in_fuel)/(O2_demand_fuel + O2_demand_third) - 1

        fuel_burn['n_air'] = n_air
        n_fuel = fuel_burn['n_fuel']
        fix_ratios(fuel_burn)
        
        # O2 excess should be wrong, same for ratios
        return fuel_burn
    if n_out is not None:
        # missing fuel and air, but know outlet and another spec
        
        # Burn all the third stream at O2 excess=0
        # Call the solver with the first
        # Add the two air inlets and the combustion products
        third_burn = fuel_air_spec_solver(zs_fuel=zs_third, n_fuel=n_third,
                                          Vm_fuel=Vm_third, MW_fuel=MW_third,
                                          combustion_stoichiometries=combustion_stoichiometries,
                                          **common_specs)
        n_out_remaining = n_out - third_burn['n_out']
        
        fuel_burn = fuel_air_spec_solver(zs_fuel=zs_fuel, n_out=n_out_remaining,
                                             Vm_fuel=Vm_fuel, MW_fuel=MW_fuel,
                                             combustion_stoichiometries=combustion_stoichiometries,
                                             **common_specs)
        
        ans = {'n_out': n_out}
        ans['n_air'] = n_air = third_burn['n_air'] + fuel_burn['n_air']
        ans['n_fuel'] = n_fuel = fuel_burn['n_fuel']
        ans['ns_out'] = [ni+nj for ni, nj in zip(third_burn['ns_out'], fuel_burn['ns_out'])]
        ans['zs_out'] = normalize(ans['ns_out'])
        ans['frac_out_O2'] = ans['ns_out'][O2_index]/n_out
        ans['frac_out_O2_dry'] = ans['ns_out'][O2_index]/(n_out - ans['ns_out'][H2O_index])
        
        
        # O2 excess definition needs to include O2 from elsewhere!
        O2_in = n_air*zs_air[O2_index] + n_fuel*zs_fuel[O2_index] + n_third*zs_third[O2_index]
        O2_demand = O2_in - ans['ns_out'][O2_index]
        ans['O2_excess'] = O2_in/O2_demand - 1
        fix_ratios(ans)
        
        return ans
        

def combustion_spec_solver(zs_air, zs_fuel, zs_third, CASs, atomss, n_third,
                           n_fuel=None, n_air=None, n_out=None,
                           O2_excess=None, frac_out_O2=None,
                           frac_out_O2_dry=None, ratio=None,
                           Vm_air=None, Vm_fuel=None, Vm_third=None, 
                           MW_air=None, MW_fuel=None, MW_third=None,
                           ratio_basis='mass', reactivities=None,
                           combustion_stoichiometries=None):
    """Solves the system of equations describing a flow of air mixing with two
    flow of combustibles, one fixed and one potentially variable, and burning
    completely. All calculated variables are returned as a dictionary.

    The variables `Vm_air`, `Vm_fuel`, `Vm_third`, `MW_air`, `MW_fuel` and 
    `MW_third` are only
    required when an air-fuel ratio is given. Howver, the ratios cannot be
    calculated for the other solve options without them.
    
    Parameters
    ----------
    zs_air : list[float]
        Mole fractions of the air; most not contain any combustibles, [-]
    zs_fuel : list[float]
        Mole fractions of the fuel; can contain inerts and/or oxygen as well,
        [-]
    zs_third : list[float]
        Mole fractions of the fixed fuel flow; can contain inerts and/or oxygen
         as well, [-]
    CASs : list[str]
        CAS numbers of all compounds, [-]
    atomss : list[dict[float]]
        List of dictionaries of elements and their counts for all molecules in
        the mixtures, [-]
    n_third : float, optional
        Flow rate of third stream, (fixed) fuel flow rate, [mol/s]
    n_fuel : float, optional
        Flow rate of fuel, [mol/s]
    n_air : float, optional
        Flow rate of air, [mol/s]
    n_out : float, optional
        Flow rate of combustion products, remaining oxygen, and inerts, [mol/s]
    O2_excess : float, optional
        The excess oxygen coming out; (O2 in)/(O2 required) - 1, [-]
    frac_out_O2 : float, optional
        The mole fraction of oxygen out, [-]
    frac_out_O2_dry : float, optional
        The mole fraction of oxygen out on a dry basis, [-]
    ratio : float, optional
        Air-fuel ratio, in the specified `basis`, [-]
    Vm_air : float, optional
        Molar volume of air, [m^3/mol]
    Vm_fuel : float, optional
        Molar volume of fuel, [m^3/mol]
    Vm_third : float, optional
        Molar volume of second fuel stream, [m^3/mol]
    MW_air : float, optional
        Molecular weight of air, [g/mol]
    MW_fuel : float, optional
        Molecular weight of fuel, [g/mol]
    MW_third : float, optional
        Molecular weight of second fuel stream, [g/mol]
    ratio_basis : str, optional
        One of 'mass', 'mole', or 'volume', [-]
    reactivities : list[bool], optional
        Optional list which can be used to mark otherwise combustible 
        compounds as incombustible and which will leave unreacted, [-]
    combustion_stoichiometries : list[dict[str, float]]
        List of return values from `combustion_stoichiometry`, can be 
        provided if precomputed [-]
        
    Returns
    -------
    results : dict
        * n_fuel : Flow rate of fuel, [mol/s]
        * n_air : Flow rate of air, [mol/s]
        * n_out : Flow rate of combustion products, remaining oxygen, and 
          inerts, [mol/s]
        * O2_excess : The excess oxygen coming out; (O2 in)/(O2 required) - 1,
          [-]
        * frac_out_O2 : The mole fraction of oxygen out, [-]
        * frac_out_O2_dry : The mole fraction of oxygen out on a dry basis, [-]
        * mole_ratio : Air-fuel mole ratio, [-]
        * mass_ratio : Air-fuel mass ratio, [-]
        * volume_ratio : Air-fuel volume ratio, [-]
        * ns_out : Mole flow rates out, [mol/s] 
        * zs_out : Mole fractions out, [-]

    Notes
    -----
    Combustion products themselves cannot be set as unreactive.
    
    The function works so long as the flow rates, molar volumes, and molecular
    weights are in a consistent basis.
    
    Handling the case of the air feed containing combustibles is not
    implemented.
    
    Examples
    --------
    
    >>> zs_air = [0.79, 0.205, 0, 0, 0, 0.0045, 0.0005]
    >>> zs_fuel = [0.025, 0.025, 0.85, 0.07, 0.029, 0.0005, 0.0005]
    >>> zs_third = [0.1, 0.005, 0.5, 0.39, 0, 0.005, 0]
    >>> CASs = ['7727-37-9', '7782-44-7', '74-82-8', '74-84-0', '74-98-6', '7732-18-5', '124-38-9']
    >>> atomss = [{'N': 2}, {'O': 2}, {'H': 4, 'C': 1}, {'H': 6, 'C': 2}, {'H': 8, 'C': 3}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}]
    >>> combustion_stoichiometries = [combustion_stoichiometry(atoms) for atoms in atomss]
    >>> ans = combustion_spec_solver(zs_air=zs_air, zs_fuel=zs_fuel, zs_third=zs_third, CASs=CASs, atomss=atomss, n_third=1.0, n_fuel=5.0, O2_excess=0.3, Vm_air=0.02493, Vm_fuel=0.02488, Vm_third=.024, MW_air=28.79341351, MW_fuel=18.55158039, MW_third=22.0)
    >>> ans['n_air']
    80.6317073170732
    """
    TRACE_FRACTION_IN_AIR = 1e-12
    if reactivities is None:
        reactivities = [True for i in zs_air]
    combustibilities = [is_combustible(CASs[i], atomss[i], reactivities[i]) for i in range(len(CASs))]

    zs_air_comb = []
    zs_air_pure = []
    n_air_comb = 0.0
    n_air_pure = 0.0

    air_has_combustibles = False
    for combustible, zi in zip(combustibilities, zs_air):
        if combustible and zi > TRACE_FRACTION_IN_AIR:
            air_has_combustibles = True
            zs_air_comb.append(zi)
            zs_air_pure.append(0.0)
            if n_air is not None:
                n_air_comb += zi*n_air
        else:
            zs_air_pure.append(zi)
            zs_air_comb.append(0.0)
            if n_air is not None:
                n_air_pure += zi*n_air
    
    if air_has_combustibles:
        try:
            zs_air_comb = normalize(zs_air_comb)
        except ZeroDivisionError:
            pass
        zs_air_pure = normalize(zs_air_pure)
        
    if not air_has_combustibles:
        return fuel_air_third_spec_solver(zs_air=zs_air, zs_fuel=zs_fuel, 
                                          zs_third=zs_third, CASs=CASs, 
                                          atomss=atomss, n_third=n_third,
                           n_fuel=n_fuel, n_air=n_air, n_out=n_out,
                           O2_excess=O2_excess, frac_out_O2=frac_out_O2,
                           frac_out_O2_dry=frac_out_O2_dry, ratio=ratio,
                           Vm_air=Vm_air, Vm_fuel=Vm_fuel, Vm_third=Vm_third, 
                           MW_air=MW_air, MW_fuel=MW_fuel, MW_third=MW_third,
                           ratio_basis=ratio_basis, reactivities=reactivities,
                           combustion_stoichiometries=combustion_stoichiometries)
    else:
        # Can handle air flow specs easily - burn the air, and then pass in the 
        # burnt air composition. This includes n_air and n_out spec.
        
        # Can handle fuel and air/fuel ratio easily as is equivalent spec to air flow.
        # Can 
        raise NotImplementedError("Composition of air includes combustibles")