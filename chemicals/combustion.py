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
from dataclasses import dataclass

__all__ = ('CombustionData',
           'CombustionStoichiometry',
           'combustion_data',
           'combustion_stoichiometry',
           'HHV_modified_Dulong',
           'HHV_stoichiometry',
           'LHV_from_HHV')

# %% Supporting functions

# TODO: Consider moving this to elements.py
def as_atoms(formula):
    if isinstance(formula, str):
        atoms = simple_formula_parser(formula)
    elif isinstance(formula, dict):
        atoms = formula
    else:
        raise ValueError("atoms must be either a string or dictionary, "
                        f"not a '{type(formula).__name__}' object")
    return atoms


# %% Combustion functions

DULONG = 'Dulong'
STOICHIOMETRY = 'Stoichiometry'
HHV_methods = (DULONG, STOICHIOMETRY)

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
    CombustionData(stoichiometry=CombustionStoichiometry(O2=-1.5, CO2=1, Br2=0.0, I2=0.0, HCl=0, HF=0, SO2=0, N2=0.0, P4O10=0.0, H2O=2.0, Ash=0), HHV=-7.26e+05, Hf=-2.39e+05, MW=32)
    
    References
    ----------
    .. [1] Green, D. W. Distillation. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018
    
    """
    atoms = as_atoms(formula)
    stoichiometry = combustion_stoichiometry(atoms)
    if not MW:
        MW = molecular_weight(atoms)
    if not method:
        method = 'Dulong' if Hf is None else 'Stoichiometry'
    if method == DULONG:
        HHV = HHV_modified_Dulong(mass_fractions(atoms, MW), MW)
        if Hf: raise ValueError("cannot specify Hf if method is 'Dulong'")
        Hf = HHV - HHV_stoichiometry(stoichiometry, 0)
    elif method == STOICHIOMETRY:
        if Hf is None: raise ValueError("must specify Hf if method is 'Stoichiometry'")
        HHV = HHV_stoichiometry(stoichiometry, Hf)
    else:
        raise ValueError("method must be either 'Stoichiometric' or 'Dulong', "
                        f"not {repr(method)}")
    return CombustionData(stoichiometry, HHV, Hf, MW)

def combustion_stoichiometry(atoms, MW=None):
    r"""
    Return a CombustionStoichiometry object that contains the cstoichiometric
    coefficients of chemical combustion as in [1]_.
    
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

    Examples
    --------
    Methane gas burning:
    
    >>> combustion_stoichiometry({'C': 1, 'H':4})
    CombustionStoichiometry(O2=-2.0, CO2=1, Br2=0.0, I2=0.0, HCl=0, HF=0, SO2=0, N2=0.0, P4O10=0.0, H2O=2.0, Ash=0)

    References
    ----------
    .. [1] Green, D. W. Distillation. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018

    """
    combustion_atoms = {i:atoms.get(i, 0) for i in combustible_elements}
    C, H, N, O, S, Br, I, Cl, F, P = combustion_atoms.values()
    MW = MW or molecular_weight(atoms)
    Ash = MW - molecular_weight(combustion_atoms)
    return CombustionStoichiometry(
        (Cl + F)/4. + O/2. - (C + S + H/4. + 5*P/4.), # O2
        C, # CO2
        Br/2., # Br2
        I/2., # I2
        Cl, # HCl
        F, # HF
        S, # SO2
        N/2., # N2
        P/4., # P4O10
        (H - Cl - F)/2., # H2O
        Ash if Ash / MW > 0.0001 else 0 # Ash
    )

def HHV_stoichiometry(stoichiometry, Hf):
    """
    Return the higher heating value [HHV; in J/mol] based on the 
    theoretical combustion stoichiometry and the heat of formation of
    the chemical as in [1]_.
    
    Parameters
    ----------
    stoichiometry : dict[str, float] or CombustionStoichiometry
        Stoichiometric coefficients of combustion. If stoichiometry is a 
        dictionary, it should inlcude the following keys: 'H2O', 'CO2', 'SO2', 
        'Br2', 'I2', 'HCl', 'Hf' and 'P4O10'.
    Hf : float
        Heat of formation [J/mol].
    
    Notes
    -----
    The combustion reaction is based on the following equation:
    
    .. math::
        C_c H_h O_o N_n S_s Br_b I_i Cl_x F_f P_p + kO_2 -> cCO_2 + \frac{b}{2}Br_2 + \frac{i}{2}I + xHCl + fHF + sSO_2 + \frac{n}{2}N_2 + \frac{p}{4}P_4O_{10} +\frac{h + x + f}{2}H_2O
        
        k = c + s + \frac{h}{4} + \frac{5P}{4} - \frac{x + f}{4} - \frac{o}{2}
    
    The HHV is calculated as the heat of reaction.
    
    References
    ----------
    .. [1] Green, D. W. Distillation. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018
    
    """
    if isinstance(stoichiometry, CombustionStoichiometry):
        stoichiometry = stoichiometry.__dict__
    return (- 285825 * stoichiometry['H2O']
            - 393474 * stoichiometry['CO2']
            - 296800 * stoichiometry['SO2']
            + 30880 * stoichiometry['Br2']
            + 62417 * stoichiometry['I2']
            - 92173 * stoichiometry['HCl']
            - 272711 * stoichiometry['HF']
            - 3009940 * stoichiometry['P4O10'] - Hf)

def HHV_modified_Dulong(mass_fractions, MW=1.):
    r"""
    Return higher heating value [HHV; in J/mol] based on the modified 
    Dulong's equation.
    
    Parameters
    ----------
    mass_fractions : dict[str, float]
        Dictionary of atomic mass fractions [-].
    MW : float
        Molecular weight of substance [g/mol].
    
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
    
    >>> HHV_modified_Dulong({'C': 0.67, 'H': 0.05, 'S': 0.015, 'N': 1.5, 'O': 0.087, 'Ash': 0.098})
    -283.75550000000004
    
    References
    ----------
    .. [1] Green, D. W. Distillation. In Perry’s Chemical Engineers’ Handbook,
       9 ed.; McGraw-Hill Education, 2018
    
    """
    C = mass_fractions.get('C', 0)
    H = mass_fractions.get('H', 0)
    O = mass_fractions.get('O', 0)
    S = mass_fractions.get('S', 0)
    if O > 0.105:
        raise ValueError("Dulong's formula is only valid at 10 wt. % Oxygen "
                        f"or less ({O:.0%} given)")
    return - MW * (338*C  + 1428*(H - O/8)+ 95*S)    

def LHV_from_HHV(HHV, N_H2O):
    """Return the lower heating value [LHV; in J/mol] of a chemical given
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


# %% Combustion reaction data types

@dataclass(frozen=True)
class CombustionStoichiometry:
    """
    Create a CombustionStoichiometry object that contains the stoichiometry 
    coefficients of the reactants and products.
    
    Parameters
    ----------
    O2: float
        Molecular oxygen [mol].
    CO2: float
        Carbon dioxide [mol].
    Br2: float
        Molecular bromine [mol].
    I2: float
        Molecular iodine [mol].
    HCl: float
        Hydrochloric acid [mol].
    HF: float
        Hydrofloric acid [mol].
    SO2: float
        Sulfur dioxide [mol].
    N2: float
        Molecular nitrogen [mol].
    P4O10: float
        Phosphorous pentoxide [mol].
    H2O: float
        Water [mol].
    Ash: float
        Other products not accounted for in stoichiometry [mol].
    
    Attributes
    ----------
    O2: float
        Molecular oxygen [mol].
    CO2: float
        Carbon dioxide [mol].
    Br2: float
        Molecular bromine [mol].
    I2: float
        Molecular iodine [mol].
    HCl: float
        Hydrochloric acid [mol].
    HF: float
        Hydrofloric acid [mol].
    SO2: float
        Sulfur dioxide [mol].
    N2: float
        Molecular nitrogen [mol].
    P4O10: float
        Phosphorous pentoxide [mol].
    H2O: float
        Water [mol].
    Ash: float
        Other products not accounted for in stoichiometry [mol].
    
    """
    O2: float
    CO2: float
    Br2: float
    I2: float
    HCl: float
    HF: float
    SO2: float
    N2: float
    P4O10: float
    H2O: float
    Ash: float

@dataclass(frozen=True)
class CombustionData:
    """
    Return a CombustionData object that contains the stoichiometry 
    coefficients of the reactants and products, the lower and higher 
    heating values [LHV, HHV; in J/mol], the heat of formation [Hf; in J/mol],
    and the molecular weight [MW; in g/mol].
    
    Parameters
    ----------
    stoichiometry : CombustionStoichiometry
        Stoichiometric coefficients of the reactants and products.
    HHV : float
        Higher heating value [J/mol].
    Hf : float
        Heat of formation [J/mol].
    MW : float
        Molecular weight [g/mol].
    
    Attributes
    ----------
    stoichiometry : CombustionStoichiometry
        Stoichiometric coefficients of the reactants and products.
    HHV : float
        Higher heating value [J/mol].
    Hf : float
        Heat of formation [J/mol].
    MW : float
        Molecular weight [g/mol].
    
    """
    stoichiometry: CombustionStoichiometry 
    HHV: float
    Hf: float
    MW: float
    
    @property
    def LHV(self):
        """Lower heating value [J/mol]."""
        return LHV_from_HHV(self.HHV, self.stoichiometry.H2O)


