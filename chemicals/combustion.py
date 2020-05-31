# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
All data and methods related to chemical combustion.
"""
from .elements import mass_fractions, molecular_weight, simple_formula_parser
from dataclasses import dataclass

__all__ = ('CombustionData',
           'CombustionStoichiometry',
           'HHV_modified_Dulong',
           'HHV_stoichiometry',
           'LHV_from_HHV')

# %% Combustion functions

combustible_elements = ('C', 'H', 'N', 'O', 'S', 'Br', 'I', 'Cl', 'F', 'P')

def combustion_data(formula, Hf=None, MW=None, method=None):
    r"""
    Return a CombustionData object that contains the stoichiometry 
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
    method : "Stoichiometry" or "Dulong"
        Method to estimate LHV and HHV.    
    
    Returns
    -------
    combustion_data : CombustionData
        Data related to combustion [-].
    
    Notes
    -----
    Default heats of formation for chemicals are at 298 K, 1 atm. The
    combustion reaction is based on the following equation:
    
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

    >>> from chemicals.combustion import CombustionData
    >>> combustion_data({'H': 4, 'C': 1, 'O': 1}, Hf=-239100)
    
    """
    
    stoichiometry = combustion_stoichiometry(formula)
    atoms = as_atoms(formula)
    if not MW:
        MW = molecular_weight(atoms)
    if not method:
        method = 'Dulong' if Hf is None else 'Stoichiometry'
    if method == DULONG:
        HHV = HHV_modified_Dulong(atoms, MW)
        if Hf: raise ValueError("cannot specify Hf if method is 'Dulong'")
        Hf = HHV - HHV_stoichiometry(stoichiometry, 0)
    elif method == STOICHIOMETRY:
        if Hf is None: raise ValueError("must specify Hf if method is 'Stoichiometry'")
        HHV = HHV_stoichiometry(stoichiometry, Hf)
    else:
        raise ValueError("method must be either 'Stoichiometric' or 'Dulong', "
                        f"not {repr(method)}")
    LHV = LHV_from_HHV(HHV, stoichiometry.H2O)
    return CombustionData(stoichiometry, HHV, LHV, Hf, MW)

def combustion_stoichiometry(atoms, MW=None):
    r"""
    Return a dictionary of the combustion stoichiometry of a chemical.
    
    Parameters
    ----------
    atoms : dict
        Dictionary of atoms and their counts.
    MW : float, optional
        Molecular weight of chemical.

    Notes
    -----
    The stoichiometry is given by:
        
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O

        k = c + s + \frac{h}{4} + \frac{5P}{4} - \frac{x + f}{4} - \frac{o}{2}

    """
    combustion_atoms = {i:atoms.get(i, 0) for i in combustible_elements}
    C, H, N, O, S, Br, I, Cl, F, P = combustion_atoms.values()
    MW = MW or molecular_weight(atoms)
    Ash = MW - molecular_weight(combustion_atoms)
    return {
        'O2': (Cl + F)/4. + O/2. - (C + S + H/4. + 5*P/4.),
        'CO2': C,
        'Br2': Br/2.,
        'I2': I/2.,
        'HCl': Cl,
        'HF': F,
        'SO2': S,
        'N2': N/2.,
        'P4O10': P/4.,
        'H2O': (H - Cl - F)/2.,
        'Ash': Ash if Ash > 0.01 else 0
    }

def HHV_stoichiometry(stoichiometry, Hf):
    """
    Return the higher heating value [HHV; in J/mol] based on the 
    theoretical combustion stoichiometry and the heat of formation of
    the chemical.
    
    Parameters
    ----------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of combustion [-]. Must inlcude the 
        following keys: 'H2O', 'CO2', 'SO2', 'Br2', 'I2', 'HCl', 'Hf' and 'P4O10'.
    Hf : float
        Heat of formation of combustable [J/mol].
    
    """
    return (- 285825 * stoichiometry['H2O']
            - 393474 * stoichiometry['CO2']
            - 296800 * stoichiometry['SO2']
            + 30880 * stoichiometry['Br2']
            + 62417 * stoichiometry['I2']
            - 92173 * stoichiometry['HCl']
            - 272711 * stoichiometry['Hf']
            - 3009940 * stoichiometry['P4O10'] - Hf)

def HHV_modified_Dulong(self, atoms, MW=None):
    r"""
    Return higher heating value [HHV; in J/mol] based on the modified 
    Dulong's equation.
    
    Parameters
    ----------
    atoms : dict[str, float]
        Dictionary of atoms and their counts.
    MW : float, optional
        Molecular weight of chemical.
    
    Notes
    -----
    The heat of combustion in J/mol is given by Dulong's equation [1]_:
    
    .. math:: 
        Hc (J/mol) = MW \cdot (338C + 1428(H - O/8)+ 95S)
    
    This equation is only good for <10 wt. % Oxygen content. Variables C, H, O,
    and S are atom weight fractions.
    
    References
    ----------
    .. [1] Brown et al., Energy Fuels 2010, 24 (6), 3639–3646.
    
    """
    massfrac = mass_fractions(atoms, MW)
    C = massfrac.get('C', 0)
    H = massfrac.get('H', 0)
    O = massfrac.get('O', 0)
    S = massfrac.get('S', 0)
    if O > 0.105:
        raise ValueError("Dulong's formula is only valid at 10 wt. % Oxygen "
                        f"or less ({O:.0%} given)")
    return - MW * (338*C  + 1428*(H - O/8)+ 95*S)    

def LHV_from_HHV(HHV, N_H2O):
    """Estimate the lower heating value [LHV; in J/mol] of a chemical given
    the higher heating value [HHV; in J/mol] and the number of water
    molecules formed per molecule burned."""
    return HHV + 44011.496 * N_H2O


# %% Combustion reaction

DULONG = 'Dulong'
STOICHIOMETRY = 'Stoichiometry'
HHV_methods = (DULONG, STOICHIOMETRY)

def as_atoms(formula):
    if isinstance(formula, str):
        atoms = simple_formula_parser(formula)
    elif isinstance(formula, dict):
        atoms = formula
    else:
        raise ValueError("atoms must be either a string or dictionary, "
                        f"not a '{type(formula).__name__}' object")
    return atoms

@dataclass(frozen=True)
class CombustionStoichiometry:
    O2: float = 0
    CO2: float = 0
    Br2: float = 0
    I2: float = 0
    HCl: float = 0
    HF: float = 0
    SO2: float = 0
    N2: float = 0
    P4O10: float = 0
    H2O: float = 0
    Ash: float = 0

@dataclass(frozen=True)
class CombustionData:
    """
    Create a CombustionData object that contains the stoichiometry 
    coefficients of the reactants and products and the lower and higher 
    heating values of a chemical [LHV, HHV; in J/mol].
    
    Parameters
    ----------
    stoichiometry : CombustionStoichiometry
        Stoichiometry coefficients of the reactants and products.
    LHV : float 
        Lower heating value [J/mol].
    HHV : float
        Higher heating value [J/mol].
    Hf : float
        Heat of formation [J/mol].
    MW : float
        Molecular weight [g/mol].
    
    """
    stoichiometry: CombustionStoichiometry 
    HHV: float
    LHV: float
    Hf: float
    MW: float
    
    def __repr__(self):
        return f"{type(self).__name__}(stoichiometry={repr(self.stoichiometry)}, HHV={self.HHV:.3g}, Hf={self.Hf:.3g}, MW={self.MW:.3g})"
    


