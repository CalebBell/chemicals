import os
from thermo import VaporPressure, Chemical
from chemicals import Tc, Pc, omega_definition
from chemicals.identifiers import pubchem_db
from joblib import Parallel, delayed

folder = os.path.join(os.path.dirname(__file__), '..', '..', 'chemicals', 'Critical Properties')

pubchem_db.autoload_main_db()
f = open(os.path.join(folder, 'omega_Psat_Tc_predictions.tsv'), 'w')

dump_oder = ['omega',]
keys = ['CAS', 'omega']
lines = ['\t'.join(keys) + '\n']
def generate_line(CASi):
    chem_info = pubchem_db.CAS_index[CASi]
    CAS = chem_info.CASs
    Tc_val = Tc(CAS)
    Pc_val = Pc(CAS)
    if Tc_val is None or Pc_val is None:
        return False
    c = Chemical(CAS)
    Psat_07 = c.VaporPressure(0.7*Tc_val)
    if Psat_07 is None or Psat_07 >= Pc_val or Psat_07 <= 0.0:
        return False
    #print(Tc_val, Psat_07, Pc_val)
    omega = omega_definition(Psat_07, Pc_val)
    
    values = [omega]
    line = values
    line.insert(0, CAS)
    for i, v in enumerate(line):
        if v is None:
            line[i] = ''
        elif v == '':
            pass
        elif isinstance(v, (int, float)):
            line[i] = '{:.8g}'.format(v)
    to_write = '\t'.join(line) + '\n'
    return to_write
    
    
#retVals = []
#for CASi in sorted(pubchem_db.CAS_index):
#    retVals.append(generate_line(CASi))

retVals = Parallel()(delayed(generate_line)(CASi) for CASi in sorted(pubchem_db.CAS_index))


for l in retVals:
    if l is not False:
        lines.append(l)

f.writelines(lines)
f.close()
