import lxml.html
import os, json, re
import appdirs
from joblib import Parallel, delayed

from fluids.constants import torr
from fluids.core import C2K, F2K, R2K
from fluids.numerics import mean

from chemicals.identifiers import serialize_formula
from chemicals.elements import molecular_weight, nested_formula_parser

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.inchi import MolToInchi, MolFromInchi, InchiToInchiKey

folder = os.path.dirname(__file__)

cache_dir = appdirs.user_cache_dir(appname='chemical_metadata')
if not os.path.exists(cache_dir):
    os.mkdir(cache_dir)

common_chemistry_cache_dir = os.path.join(cache_dir, 'common_chemistry')
if not os.path.exists(common_chemistry_cache_dir):
    os.mkdir(common_chemistry_cache_dir)

CASs_process = open(os.path.join(folder, 'common_chemistry_CASs.tsv')).readlines()
CASs_process = [i.strip() for i in CASs_process]


def remove_html(text):
    '''
    Examples
    --------
    
    >>> remove_html('Boric acid (H<sub>3</sub>BO<sub>3</sub>), trimethyl ester')
    'Boric acid (H3BO3), trimethyl ester'
    
    '''
    if '<' in text:
        return lxml.html.fromstring(text).text_content()
    return text

number_range = re.compile(r'([+-]?(\d*\.)?\d+)-([+-]?(\d*\.)?\d+)')

def temperature_if_present(s):
    '''Identify a temperature if one is in the
    string provided, and return it in Kelvin units.
    If a range is given, this returns the average.
    '''
    if '@ Temp' not in s or '°C' not in s:
        return None
    number_string = s.split('@ Temp')[1].split('°C')[0]
    if re.search(number_range, number_string):
        # If there are two numbers, split on the range
        bits = number_range.split(number_string)
        low, high = bits[1], bits[3]
        temp = mean([float(v) for v in (low, high)])
    else:
        temp = float(number_string.replace(':', '').strip())
    temp = C2K(temp)
    return temp


def pressure_if_present(s):
    '''Identify a pressure if one is in the
    string provided, and return it in Pa units.
    If a range is given, this returns the average.
    '''
    if 'Press' not in s or 'Torr' not in s:
        return None
    number_string = s.split(' Press')[1].split('Torr')[0]
    if re.search(number_range, number_string):
        # If there are two numbers, split on the range
        bits = number_range.split(number_string)
        low, high = bits[1], bits[3]
        pres = mean([float(v) for v in (low, high)])
    else:
        pres = float(number_string.replace(':', '').strip())
    pres = float(pres)*torr
    return pres

def process_temperataure(s):
    '''Process a string containing either a 
    melting point or a boiling point. These often have
    a pressure attatched.
    '''
    if '°C' not in s:
        print(s)
        raise ValueError('Unknown units')
    number = s.split('°C')[0]
    P = pressure_if_present(s)
    if re.search(number_range, number):
        # If there are two numbers, split on the range
        bits = number_range.split(s)
        low, high = bits[1], bits[3]
        number = mean([float(v) for v in (low, high)])
    else:
        number = float(number.replace(':', '').strip())
    number = C2K(number)
    return number, P

def process_density(s):
    '''Process a string containing a mass density.
    These often have a pressure and/or temperature attatched.
    '''
#     print(s, 'in density string')
    if 'g/cm<sup>3</sup' not in s and 'g/cm3' not in s:
        print(s)
        raise ValueError
    if 'g/cm<sup>3</sup' in s:
        start, *end = s.split('g/cm<sup>3</sup>')
    else:
        start, *end = s.split('g/cm3')

    if re.search(number_range, s):
        # If there are two numbers, split on the range
        bits = number_range.split(s)
        low, high = bits[1], bits[3]
        start = mean([float(v) for v in (low, high)])
        
    temp = temperature_if_present(s)
    P = pressure_if_present(s)
    rho = float(start)*1000
#     print(rho, temp, P)
    return rho, temp, P
    
    
cases = [('250-257 °C', (526.65, None)), 
         ('-50 °C', (223.14999999999998, None)),
        ('100 °C', (373.15, None))]
for t_str, expect in cases:
    assert process_temperataure(t_str) == expect

cases = [('1.1296 g/cm<sup>3</sup>', (1129.6, None, None)),
        ('1.0818 g/cm3 @ Temp: 41 °C', (1081.8000000000002, 41+273.15, None)),
        ('2.09-2.23 g/cm3 ', (2160, None, None)),
        ('0.944-0.952 g/cm3 @ Temp: 20-25 °C', (948, 22.5+273.15, None)),
        ('0.8495 g/cm3 @ Temp: 25 °C; Press: 4 Torr', (849.5, 298.15, 4*torr))]
for r_str, expect in cases:
    assert process_density(r_str) == expect

def MW_from_smiles(smiles):
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            formula = CalcMolFormula(mol, True, True)
            formula = serialize_formula(formula)
            MW = molecular_weight(nested_formula_parser(formula))
            return MW
    
def common_chemistry_data(CASRN):
    '''Load the chemical data for the specified CAS.
    This function returns None if no data is available.
    '''
    file_path = os.path.join(common_chemistry_cache_dir, CASRN)
    file_data = open(file_path).read()
    json_data = json.loads(file_data)
#     print(json_data)

    if 'message' in json_data and 'Detail not found' == json_data['message']:
        return None

    synonyms = [remove_html(k) for k in json_data['synonyms']]
    inchi = json_data['inchi'] if json_data['inchi']  else None
    inchiKey = json_data['inchiKey'] if json_data['inchiKey']  else None
    smiles = json_data['canonicalSmile'] if json_data['canonicalSmile']  else None
    formula = remove_html(json_data['molecularFormula']) if json_data['molecularFormula']  else None 
    MW = float(json_data['molecularMass']) if json_data['molecularMass'] else None
    #if MW is None:
        #MW = MW_from_smiles(smiles)
        
    Tm, Tb, TbP, rho, rhoT, rhoP = None, None, None, None, None, None
    
    for prop in json_data['experimentalProperties']:
        orig_prop_val = prop['property']
        prop_val = orig_prop_val # remove_html
        prop_val = prop_val.replace('g/cm<sup>3</sup>', 'g/cm3')
        prop_val = prop_val.replace('&gt;', '>')
        prop_val = prop_val.replace('&lt;', '<')
        prop_val = prop_val.replace(' x 10<sub>', 'E-').replace('</sub>', '')
        prop_val = prop_val.replace(' x 10<sup>', 'E').replace('</sup>', '')
        prop_val = prop_val.replace(' x 10', 'E') # handle scientific notation
        if '>' in prop_val or '<' in prop_val:
            # Not a real property
            continue
        if prop['name'] == 'Boiling Point':
            Tb, TbP = process_temperataure(prop_val)
        if prop['name'] == 'Melting Point':
            Tm, _ = process_temperataure(prop_val)
        if prop['name'] == 'Density':
            rho, rhoT, rhoP = process_density(prop_val)
    return {'CAS': json_data['rn'], 'name': remove_html(json_data['name']), 'synonyms': synonyms,
            'formula': formula, 'inchi': inchi, 'inchiKey': inchiKey, 'smiles': Chem,
           'Tm': Tm, 'Tb': Tb, 'TbP': TbP, 'rho': rho, 'rhoT': rhoT, 'rhoP': rhoP, 'MW': MW,
           }
'''Load and process the data into dictionaries in-memory.
'''
processed_data = Parallel(n_jobs=16, batch_size=10000, verbose=60)(delayed(common_chemistry_data)(c) for c in CASs_process)
print('Loaded %d chemicals from specified CAS list' %(len(processed_data)))

chemical_data = {}
for dat in processed_data:
    if dat['Tb'] is not None or dat['Tm'] is not None or dat['rho'] is not None:
        chemical_data[dat['CAS']] = dat
print('%d chemicals have at least one experimental property' %(len(chemical_data)))
