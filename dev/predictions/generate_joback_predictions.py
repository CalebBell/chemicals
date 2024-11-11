import os

from joblib import Parallel, delayed
from rdkit import Chem
from thermo import Joback

from chemicals.identifiers import pubchem_db

folder = os.path.join(os.path.dirname(__file__), '..', '..', 'chemicals', 'Misc')

f = open(os.path.join(folder, 'joback_predictions.tsv'), 'w')
from rdkit import Chem
pubchem_db.finish_loading()
dump_oder = ['Tm', 'Hfus', 'Hvap','Tb', 'Tc', 'Pc', 'Vc', 'Hf']
keys = ['CAS', 'Tm', 'Hfus', 'Hvap','Tb', 'Tc', 'Pc', 'Vc', 'Hfg', 'Cpg0', 'Cpg1', 'Cpg2', 'Cpg3', 'mul0', 'mul1']
lines = ['\t'.join(keys) + '\n']
def generate_line(CASi):
    pubchem_db.finish_loading() # needed for multiprocessing?
    chem_info = pubchem_db.CAS_index[CASi]
    mol = Chem.MolFromSmiles(chem_info.smiles)
    if mol is None:
        # rdkit coult not parse
        return False
    obj = Joback(mol)
    if obj.success:
        try:
            estimates = obj.estimate()
        except:
            # Did not match any groups
            return False
        values = [estimates[k] for k in dump_oder]
        line = values
        line.insert(0, CASi)
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
            elif i == 0:  # CAS number - keep as is
                line[i] = str(v)
            elif v == '':
                pass
            elif isinstance(v, (int, float)):
                line[i] = f'{v:.8g}'
        to_write = '\t'.join(line) + '\n'
        return to_write
    return False



retVals = Parallel(n_jobs=16, batch_size=500, backend='multiprocessing')(delayed(generate_line)(CASi) for CASi in sorted(pubchem_db.CAS_index))
for l in retVals:
    if l is not False:
        lines.append(l)

f.writelines(lines)
f.close()
