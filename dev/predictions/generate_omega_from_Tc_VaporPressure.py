import os

from joblib import Parallel, delayed
from thermo import Chemical

from chemicals import Pc, Tc, omega_definition
from chemicals.identifiers import pubchem_db

folder = os.path.join(os.path.dirname(__file__), '..', '..', 'chemicals', 'Critical Properties')

pubchem_db.autoload_main_db()

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

    # Add some basic limits for omega
    if omega > 3.0:
        omega = 3.0
    if omega < -2.0:
        omega = -2.0

    values = [omega]
    line = values
    line.insert(0, CASi)
    for i, v in enumerate(line):
        if v is None:
            line[i] = ''
        elif v == '':
            pass
        elif i == 0:
            line[i] = str(v) # CAS
        elif isinstance(v, (int, float)):
            line[i] = f'{v:.8g}'
    to_write = '\t'.join(line) + '\n'
    return to_write


#retVals = []
#for CASi in sorted(pubchem_db.CAS_index):
#    retVals.append(generate_line(CASi))

retVals = Parallel()(delayed(generate_line)(CASi) for CASi in sorted(pubchem_db.CAS_index))


for l in retVals:
    if l is not False:
        lines.append(l)
f = open(os.path.join(folder, 'omega_Psat_Tc_predictions.tsv'), 'w')
f.writelines(lines)
f.close()
