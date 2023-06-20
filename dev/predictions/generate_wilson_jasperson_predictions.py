import os

from rdkit import Chem
from thermo import Wilson_Jasperson

from chemicals import Tb
from chemicals.identifiers import pubchem_db

folder = os.path.join(os.path.dirname(__file__), '..', '..', 'chemicals', 'Critical Properties')
nan = float("nan")
pubchem_db.autoload_main_db()
f = open(os.path.join(folder, 'wilson_jasperson_Tc_Pc_predictions.tsv'), 'w')
from rdkit import Chem

dump_oder = ['Tc', 'Pc']
keys = ['CAS', 'Tc', 'Pc']
lines = ['\t'.join(keys) + '\n']
def generate_line(CASi):
    chem_info = pubchem_db.CAS_index[CASi]
    CAS = chem_info.CASs
    mol = Chem.MolFromSmiles(chem_info.smiles)
    if mol is None:
        # rdkit coult not parse
        return False
    Tb_db = Tb(CAS)
    if Tb_db is None:
        return False
    Tc, Pc, missing_Tc_increments, missing_Pc_increments = Wilson_Jasperson(mol, Tb=Tb_db)
    if missing_Tc_increments:
        Tc = nan
    if missing_Pc_increments:
        Pc = nan
    if not missing_Tc_increments or not missing_Pc_increments:
        values = [Tc, Pc]
        line = values
        line.insert(0, CASi)
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


for CASi in sorted(pubchem_db.CAS_index):
    ans = generate_line(CASi)
    if ans is not False:
        lines.append(ans)

#retVals = Parallel()(delayed(generate_line)(CASi) for CASi in sorted(pubchem_db.CAS_index))
#for l in retVals:
#    if l is not False:
#        lines.append(l)

f.writelines(lines)
f.close()
