"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2022, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
SOFTWARE.
"""

import os

import numpy as np
import pandas as pd
from thermo import *

from chemicals import *
from chemicals import data_reader
from chemicals.identifiers import pubchem_db
from chemicals.utils import source_path

pubchem_db.finish_loading()

# Can't use the database
data_reader.USE_CONSTANTS_DATABASE = False
_RI = RI
def RI(CASRN):
    return _RI(CAS)[0]

def RIT(CASRN):
    return _RI(CAS)[1]

funcs = [MW, Tt, Tm, Tb, Tc, Pt, Pc, Vc,
         Zc,omega, T_flash, T_autoignition, LFL, UFL,
                     Hfs, Hfl, Hfg, S0s, S0l, S0g,
                     Hfus,  Stockmayer, molecular_diameter,
                     dipole_moment, logP, RG, RON, MON, ignition_delay, linear,
                     GWP, ODP, RI, RIT, GTP, 
                     hansen_delta_d, hansen_delta_p, hansen_delta_h]
props = [f.__name__ for f in funcs]
CASs = []
prop_array = [[] for _ in range(len(props))]

objs = list(pubchem_db.CAS_index.values())
stuff = []
hits = 0
for v in objs:
    CAS = v.CASs
    CASs.append(CAS)
    for i, p in enumerate(props):
        if p == 'MW':
            lookedup_constant = float(v.MW)
        else:
            lookedup_constant = funcs[i](CASRN=CAS)
        if lookedup_constant is None:
            lookedup_constant = float("nan")
        prop_array[i].append(lookedup_constant)



prop_array_T = np.array(prop_array).T
CASs_ints = [CAS_to_int(i) for i in CASs]

# Would not be good if there were multiple values
assert len(CASs_ints) == len(set(CASs_ints))


df = pd.DataFrame(prop_array_T, columns=props, index=CASs_ints)

# Does not save memory except when compressed and causes many test failures
#df = df.fillna(value=np.nan).astype(np.float32)
# df.sort_values(by=['logP'], inplace=True)
df.sort_index(inplace=True)

from sqlalchemy import create_engine

print('going to dump to database now')
db_dest = os.path.join(source_path, 'Misc', 'default.sqlite')
if os.path.exists(db_dest):
    os.remove(db_dest)
engine = create_engine('sqlite:///' + db_dest, echo=False)
print('calling pandas to_sql')
df.to_sql('constants', con=engine)
