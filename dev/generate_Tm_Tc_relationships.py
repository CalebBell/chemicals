"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2021, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

import numpy as np
from thermo import *

from chemicals import *
from chemicals.critical import critical_data_CRC, critical_data_IUPAC, critical_data_Matthews, critical_data_PassutDanner

dfs = [critical_data_IUPAC, critical_data_Matthews, critical_data_CRC, critical_data_PassutDanner]
inorganic_Tc_CASs = []
for d in dfs:
    index = d.index
    for CAS in index:
        c = None
        try:
            c = Chemical(CAS)
        except:
            continue
        mol = c.rdkitmol_Hs
        if mol is not None:
            inorganic_Tc_CASs.append(CAS)

inorganic_Tc_CASs = list(set(inorganic_Tc_CASs))
len(inorganic_Tc_CASs)


filters = {'organic': lambda c: is_organic(c.rdkitmol_Hs),
'inorganic': lambda c: is_inorganic(c.rdkitmol_Hs),
'element': lambda c:  is_inorganic(c.rdkitmol_Hs) and len(c.atoms) == 1,
'binary': lambda c:  is_inorganic(c.rdkitmol_Hs) and len(c.atoms) == 2,
'ternary': lambda c:  is_inorganic(c.rdkitmol_Hs) and len(c.atoms) == 3,
'Cl': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'Cl' in c.atoms,
'N': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'N' in c.atoms,
'C': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'C' in c.atoms,
'F': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'F' in c.atoms,
'O': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'O' in c.atoms,
'I': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'I' in c.atoms,
'Br': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'Br' in c.atoms,
'Si': lambda c: is_inorganic(c.rdkitmol_Hs) and len(c.atoms) >1 and 'Si' in c.atoms,

}

for method, lbd in filters.items():
    working_Tcs, working_Tms, working_CASs, working_chemicals, working_MWs = [], [], [], [], []
    for CAS in inorganic_Tc_CASs:
        crit = Tc(CAS)
        boil = Tm(CAS)
        c = Chemical(CAS)
        if crit is not None and boil is not None and lbd(c):
            working_Tcs.append(crit)
            working_Tms.append(boil)
            working_CASs.append(CAS)
            working_chemicals.append(c)
            working_MWs.append(c.similarity_variable)

    Tc_ratios = np.array(working_Tcs)/np.array(working_Tms)
    print(method)
    print(len(working_chemicals), round(np.mean(Tc_ratios),3), np.min(Tc_ratios), np.max(Tc_ratios))
