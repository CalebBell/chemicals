import os, json
import copy
import pandas as pd
from chemicals.utils import to_num, Vm_to_rho, rho_to_Vm

emptydict = {"Name":None, "MW":None, "Tc":None, "T":[], "P":[], "Density (l)":[],
             "Density (g)":[], "Hvap":[], "Cp (l)":[], "Cp (g)":[],
            "Mu (l)":[], "Mu (g)":[], "K (l)":[], "K (g)":[], "Pr (l)":[],
            "Pr (g)":[], "sigma":[], "Beta":[],
            "Volume (l)":[], "Volume (g)":[]}

_VDISaturationDict = {}
with open('VDI Saturation Compounds Data.csv') as f:
    next(f)
    for line in f:
        values = to_num(line.strip('\n').split('\t'))
        (CASRN, _name, _MW, _Tc, T, P, rhol, rhog, Hvap, cpl, cpg, mul, mug, kl, kg, prl, prg, sigma, Beta) = values
        newdict = (_VDISaturationDict[CASRN] if CASRN in _VDISaturationDict else copy.deepcopy(emptydict))
        newdict["Name"] = _name
        newdict["MW"] = _MW
        newdict["Tc"] = _Tc
        newdict["T"].append(T)
        newdict["P"].append(P)
        newdict["Density (l)"].append(rhol) # Not actually used
        newdict["Density (g)"].append(rhog) # Not actually used
        newdict["Hvap"].append(Hvap)
        newdict["Cp (l)"].append(cpl) # Molar
        newdict["Cp (g)"].append(cpg) # Molar
        newdict["Mu (l)"].append(mul)
        newdict["Mu (g)"].append(mug)
        newdict["K (l)"].append(kl)
        newdict["K (g)"].append(kg)
        newdict["Pr (l)"].append(prl)
        newdict["Pr (g)"].append(prl)
        newdict["sigma"].append(sigma)
        newdict["Beta"].append(Beta)
        newdict["Volume (l)"].append(rho_to_Vm(rhol, _MW))
        newdict["Volume (g)"].append(rho_to_Vm(rhog, _MW))
        _VDISaturationDict[CASRN] = newdict

current_folder = os.path.dirname(__file__)
print(current_folder)
dest_pth = os.path.join(current_folder, '..', '..', 'chemicals/Misc/VDI Saturation Compounds Data.json')
print(dest_pth)
f = open(dest_pth, 'w')
json.dump(_VDISaturationDict, f, indent=2,separators=(',', ':'))
f.close()
