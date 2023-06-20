import requests
from wikidata.client import Client
from wikidata.datavalue import DatavalueError

client = Client()

import os
from time import time

from fluids import C2K, F2K, R2K
from fluids.numerics import mean
from joblib import Parallel, delayed

from chemicals import molecular_weight, nested_formula_parser
from chemicals.identifiers import CAS_to_int
from chemicals.safety import mgm3_to_ppmv

"""These are the propertiese and string representations we will retrieve"""
properties = {'P2101': 'Tm',
              'P2102': 'Tb',
              'P1109': 'RI',
              'P2066': 'Hfus',
              'P2116': 'Hvap',
              #'P2054': 'rho',
              'P2202': 'LFL',
              'P2203': 'UFL',
              'P2128': 'T_flash',
              'P2199': 'T_autoignition',
              'P3078': 'Hf', # Standard state only
              'P3071': 'S0',
              'P2993': 'logP',
             }

max_chemicals = 100000 # For testing; set to a large number to ignore this; set to maybe 100 to test a small query.
# max_chemicals = 30

wiki_prop_IDs = list(properties.keys())
wiki_prop_items = [client.get(v, load=True) for v in wiki_prop_IDs]

class NoUnitsError(Exception):
    """Generic error to raise when not enough values are given."""

CAS_item = client.get('P231', load=True)
formula_item = client.get('P274', load=True)
smiles_item = client.get('P233', load=True)
inchi_item = client.get('P234', load=True)
inchikey_item = client.get('P235', load=True)
MW_item = client.get('P2067', load=True)

degC_item = client.get('Q25267', load=True)
degF_item = client.get('Q42289', load=True)
degR_item = client.get('Q37732658', load=True)
vol_percent_item = client.get('Q2080811', load=True)
K_item = client.get('Q11579', load=True)
kJ_mol_item = client.get('Q752197', load=True)
J_mol_item = client.get('Q13035094', load=True)
J_mol_K_item = client.get('Q20966455', load=True)
g_m3_item = client.get('Q21604951', load=True)

def value_to_SI(value, unit, value_key, MW):
    if value_key in ('Tm', 'Tb', 'T_flash', 'T_autoignition'):
        if unit is degC_item:
            return C2K(value)
        if unit is degF_item:
            return F2K(value)
        if unit is degR_item:
            return R2K(value)
        if unit is K_item:
            return value
    if value_key in ('UFL', 'LFL'):
        if unit is vol_percent_item:
            return value* .01
        if unit is g_m3_item:
            return 1e-6*mgm3_to_ppmv(value, MW)
    if value_key in ('RI', 'logP'):
        return value
    if value_key in ('Hvap', 'Hfus', 'Hf'):
        if unit is kJ_mol_item:
            return value * 1e3
        if unit is J_mol_item:
            return value
    if value_key in ('S0',):
        if unit is J_mol_K_item:
            return value
    """Need to add more units or fix an item in wikidata when this happens"""
    print(value, unit, value_key)
    raise NoUnitsError('No Units')

subscripts = '₀₁₂₃₄₅₆₇₈₉'
numbers = '0123456789'
translate_subscripts = str.maketrans(subscripts, numbers)



chemical_IDs = set()
start = time()
# There isn't a need to restrict to our search to chemicals specifically; anything with a CAS will do. Some chemicals aren't mapped as chemicals.
"""
        {?cmpnd wdt:P279 wd:Q11173 .} UNION
        {?cmpnd wdt:P31 wd:Q11173 .} UNION
        {?cmpnd wdt:P31 wd:Q11344 .}

"""
wikidata_url = 'https://query.wikidata.org/sparql'
query = """
SELECT * WHERE {
  {?cmpnd wdt:P231 ?cas .}
"""
prop_keys = list(properties.keys())
for i in range(len(prop_keys)-1):
    query += '  {?cmpnd wdt:%s ?property_vaue .} UNION\n' %(prop_keys[i])
query += '  {?cmpnd wdt:%s ?property_vaue .}\n' %(prop_keys[-1])

query += '}'
r = requests.get(wikidata_url, params={'format': 'json', 'query': query})
data = r.json()
found_IDs = []
for d in data['results']['bindings']:
    s = d['cmpnd']['value']
    found_IDs.append(s.replace('http://www.wikidata.org/entity/', ''))
chemical_IDs.update(found_IDs)
# print(chemical_IDs)
# print(list(chemical_IDs)[-1], type(list(chemical_IDs)[-1]))
print(f'Identified {len(chemical_IDs)} chemicals in {time() - start} seconds')# 9344





start = time()

def process_single(quantity, value_key, MW):
    value = quantity.amount
    value = value_to_SI(value, quantity.unit, value_key, MW)
    return value

# Silver, iron, copper, arsenic, gold all have time decode errors
bad_IDs = {'Q1090', 'Q677', 'Q753', 'Q871', 'Q897'}
bad_IDs = set()
chemical_data = {}
print('Processing %d Chemicals' %(len(chemical_IDs)))

processed_items = 0
def retrieve_chemical_data(ID):
    global processed_items

    entity = client.get(ID, load=True)

    name = str(entity.label)
    CAS = entity[CAS_item]
    formula = entity.get(formula_item, None)
    if formula is not None:
        formula = formula.translate(translate_subscripts)
    smiles = entity.get(smiles_item, None)
    inchi = entity.get(inchi_item, None)
    inchikey = entity.get(inchikey_item, None)
    MW = entity.get(MW_item, None)
    if MW is not None:
        MW = MW.amount
    elif type(formula) is str:
        try:
            MW = molecular_weight(nested_formula_parser(formula))
        except:
            MW = None
    else:
        MW = None
    assert  type(CAS_to_int(CAS)) is int
    # print('CAs', type(CAS), CAS_to_int(CAS), type(CAS_to_int(CAS)), 'hi')
    dat = {'CAS': CAS_to_int(CAS), 'name': name, 'MW': MW,
                          'formula': formula, 'smiles': smiles,
                          'inchi': inchi, 'inchikey': inchikey}

    for tn, ti in zip(wiki_prop_IDs, wiki_prop_items):
        prop_key = properties[tn]
        value = None
        found = False
        if ti in entity:
            try:
                value = process_single(entity[ti], prop_key, MW)
            except NoUnitsError:
                pass
#             print(entity, value)
            if value is not None:
                dat[prop_key] = value
                found = True
        if not found:
            try:
                for item_type, value_list in entity.iterlists():
                    if item_type == ti:
                        vals = []
                        for thing in value_list:
                            try:
                                v = process_single(thing, prop_key, MW)
                                vals.append(v)
                            except NoUnitsError:
                                pass

                        if len(vals):
                            value = mean(vals)
                            dat[prop_key] = value
            except DatavalueError:
                print(ID)
                continue
    processed_items += 1
    if (processed_items %100) == 0:
        print('Processed %d items' %(processed_items))
    return dat

res =Parallel(n_jobs=16, prefer="threads")(delayed(retrieve_chemical_data)(c) for c in list(sorted(chemical_IDs)[0:max_chemicals]))
for d in res:
    chemical_data[d['CAS']] = d

print(f"Processed {len(chemical_data)} items in {time() - start} seconds")




#keys = ['CAS', 'name', 'formula', 'MW', 'Tm', 'Tb', 'Hfus', 'Hvap', 'Hf', 'S0', 'LFL', 'UFL', 'T_flash', 'T_autoignition', 'RI', 'RIT', 'logP']
keys = ['CAS', 'Tm', 'Tb', 'Hfus', 'Hvap', 'Hf', 'S0', 'LFL', 'UFL', 'T_flash', 'T_autoignition', 'RI', 'RIT', 'logP']
lines = []
for CAS in sorted(chemical_data.keys()):
    d = chemical_data[CAS]
    values = [d.get(k, '') for k in keys]
    for i, v in enumerate(values):
        if i == 0:
            # CAS`
            values[i] = str(v)
        elif v is None:
            values[i] = ''
        elif v == '':
            pass
        elif isinstance(v, (int, float)):
            values[i] = f'{v:.8g}'
        elif type(v) is str:
            pass
        else:
            values[i] = str(v)
        v = values[i]
        v = v.replace('\t', ' ')
    to_write = '\t'.join(values) + '\n'
    lines.append(to_write)

lines.sort() # history
lines.insert(0, '\t'.join(keys) + '\n')

folder = os.path.dirname(__file__)
f = open(os.path.join(folder, '..', '..', 'chemicals', 'Misc', 'wikidata_properties.tsv'), 'w')
f.writelines(lines)
f.close()
