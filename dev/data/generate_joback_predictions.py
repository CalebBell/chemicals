import os
from thermo.joback import *
from chemicals.identifiers import pubchem_db
import rdkit
from rdkit import Chem

folder = os.path.join(os.path.dirname(__file__), '..', '..', 'chemicals', 'Misc')

pubchem_db.autoload_main_db()
f = open(os.path.join(folder, 'joback_predictions.tsv'), 'w')
from rdkit import Chem

dump_oder = ['Tm', 'Hfus', 'Hvap','Tb', 'Tc', 'Pc', 'Vc', 'Hf']
keys = ['CAS', 'Tm', 'Hfus', 'Hvap','Tb', 'Tc', 'Pc', 'Vc', 'Hfg', 'Cpg0', 'Cpg1', 'Cpg2', 'Cpg3', 'mul0', 'mul1']
lines = ['\t'.join(keys) + '\n']

for CASi in sorted(pubchem_db.CAS_index):
    chem_info = pubchem_db.CAS_index[CASi]
    CAS = chem_info.CASs
    mol = Chem.MolFromSmiles(chem_info.smiles)
    if mol is None:
        # rdkit coult not parse
        continue
    obj = Joback(mol)
    if obj.success:
        try:
            estimates = obj.estimate()
        except:
            # Did not match any groups
            continue
        values = [estimates[k] for k in dump_oder]
        line = values
        line.insert(0, CAS)
        if estimates['Cpig_coeffs'] is not None:
            line.extend(estimates['Cpig_coeffs'])
        else:
            line.extend(['', '', '', ''])
        if estimates['mul_coeffs'] is not None:
            line.extend(estimates['mul_coeffs'])
        else:
            line.extend(['', ''])
        for i, v in enumerate(line):
            if v is None:
                line[i] = ''
            elif v == '':
                pass
            elif isinstance(v, (int, float)):
                line[i] = '{:.8g}'.format(v)
        to_write = '\t'.join(line) + '\n'
        lines.append(to_write)
f.writelines(lines)
f.close()
