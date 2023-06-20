import os

from joblib import Parallel, delayed
from rdkit import Chem
from thermo import Fedors

from chemicals.identifiers import pubchem_db

folder = os.path.join(os.path.dirname(__file__), '..', '..', 'chemicals', 'Critical Properties')

pubchem_db.autoload_main_db()
f = open(os.path.join(folder, 'fedors_Vc_predictions.tsv'), 'w')
from rdkit import Chem

dump_oder = ['Vc',]
keys = ['CAS', 'Vc']
lines = ['\t'.join(keys) + '\n']
def generate_line(CASi):
    chem_info = pubchem_db.CAS_index[CASi]
    CAS = CASi
    mol = Chem.MolFromSmiles(chem_info.smiles)
    if mol is None:
        # rdkit coult not parse
        return False
    Vc, status, _, _, _ = Fedors(mol)
    if status == 'OK':
        values = [Vc]
        line = values
        line.insert(0, CAS)
        for i, v in enumerate(line):
            if v is None:
                line[i] = ''
            elif v == '':
                pass
            elif isinstance(v, (int, float)):
                line[i] = f'{v:.8g}'
        to_write = '\t'.join(line) + '\n'
        return to_write
    return False



retVals = Parallel()(delayed(generate_line)(CASi) for CASi in sorted(pubchem_db.CAS_index))
for l in retVals:
    if l is not False:
        lines.append(l)

f.writelines(lines)
f.close()
