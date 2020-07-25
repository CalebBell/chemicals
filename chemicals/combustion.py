# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.'''
from chemicals.elements import mass_fractions, molecular_weight, simple_formula_parser

__all__ = ('combustion_stoichiometry',
           'CombustionData',
           'combustion_data',
           'HHV_modified_Dulong',
           'HHV_stoichiometry',
           'LHV_from_HHV')


# TODO: Consider moving this to elements.py
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
        reaction below. If 'elemental', returns those atoms in the monatomic 
        state; if 'Ash', converts all missing attoms to 'Ash' in the output at
        a `MW` of 1 g/mol, [-]

    Returns
    -------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of combustion. May inlcude the following 
        keys for complete combustion: 'H2O', 'CO2', 'SO2', 'Br2', 'I2', 'HCl', 
        'HF' 'P4O10'; if `missing_handling` is 'elemental' can include the 
        other elements; if `missing_handling` is 'Ash', Ash will be present in 
        the output if the compounds whose reactions are not included here.
        'O2' is always present, with negative values indicating oxygen is
        required. [-]

    Notes
    -----
    The stoichiometry is given by:
        
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O

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
    {'O2': -2.0, 'CO2': 1, 'H2O': 2.0}

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

    products['O2'] = -(nC + nS + .25*nH + 1.25*nP - .25*(nCl + nF) - .5*nO)
    nCO2 = nC
    nBr2 = .5*nBr
    nI2 = .5*nI

    nHCl = nCl
    nHF = nF

    nSO2 = nS

    nN2 = .5*nN
    nP4O10 = .25*nP
    nH2O = 0.5*(nH - nCl - nF)
    
    if nCO2 !=  0.0:
        products['CO2'] = nCO2
    if nBr2 != 0.0:
        products['Br2'] = nBr2
    if nI2 != 0.0:
        products['I2'] = nI2
    if nCl != 0.0:
        products['HCl'] = nCl
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
    C = mass_fractions.get('C', 0)
    H = mass_fractions.get('H', 0)
    O = mass_fractions.get('O', 0)
    S = mass_fractions.get('S', 0)
    if O > 0.105:
        raise ValueError("Dulong's formula is only valid at 10 wt. % Oxygen "
                         "or less (%s given)" %(O))
    return - (338*C  + 1428*(H - O/8)+ 95*S)    

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
        
        H_{vap} = 44011.496 \frac{J}{mol H_2O}
        
        H_2O = \frac{mol H_2O}{mol}
    
    Examples
    --------
    Methanol lower heat of combustion:
    
    >>> LHV_from_HHV(-726024.0, 2)
    -638001.008
    
    """
    return HHV + 44011.496 * N_H2O

def combustion_data(formula, Hf=None, MW=None, method=None):
    r"""
    Return a CombustionData object (a named tuple) that contains the stoichiometry 
    coefficients of the reactants and products, the lower and higher 
    heating values [LHV, HHV; in J/mol], the heat of formation [Hf; in J/mol],
    and the molecular weight [MW; in g/mol].

    Parameters
    ----------
    formula : str, or dict[str, float]
        Chemical formula as a string or a dictionary of atoms and their counts.
    MW : float, optional
        Molecular weight of chemical [g/mol].
    Hf : float, optional
        Heat of formation of given chemical [J/mol].
        Required if method is "Stoichiometry".
    method : "Stoichiometry" or "Dulong", optional
        Method to estimate LHV and HHV.    
    
    Returns
    -------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of combustion. May inlcude the following 
        keys: 'H2O', 'CO2', 'SO2', 'Br2', 'I2', 'HCl', 'HF' and 'P4O10'.
    HHV : float
        Higher heating value [J/mol].
    Hf : float
        Heat of formation [J/mol].
    MW : float
        Molecular weight [g/mol].
    
    Notes
    -----
    The combustion reaction is based on the following equation:
    
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O
        
        k = c + s + \frac{h}{4} + \frac{5P}{4} - \frac{x + f}{4} - \frac{o}{2}

    If the method is "Stoichiometry", the HHV is found using 
    through an energy balance on the reaction (i.e. heat of reaction).
    If the method is "Dulong", Dulong's equation is used [1]_:

    .. math:: 
        Hc (J/mol) = MW \cdot (338C + 1428(H - O/8)+ 95S)
        
    The LHV is calculated as follows:
        
    .. math::
        LHV = HHV + H_{vap} \cdot H_2O
        
        H_{vap} = 44011.496 \frac{J}{mol H_2O}
        
        H_2O = \frac{mol H_2O}{mol}
        
    Examples
    --------
    Liquid methanol burning:

    >>> combustion_data({'H': 4, 'C': 1, 'O': 1}, Hf=-239100)
    CombustionData(stoichiometry={'O2': -1.5, 'CO2': 1, 'H2O': 2.0}, HHV=-726024.0, Hf=-239100, MW=32.04186)
    
    References
    ----------
    .. [1] Green, D. W. Waste management. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018
    
    """
    atoms = as_atoms(formula)
    stoichiometry = combustion_stoichiometry(atoms)
    if not MW:
        MW = molecular_weight(atoms)
    if not method:
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
        return LHV_from_HHV(self.HHV, self.stoichiometry['H2O'])
