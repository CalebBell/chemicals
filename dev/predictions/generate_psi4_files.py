'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2022 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
'''

import json
import os

from fluids.constants import c, h, pi
from thermo import *

from chemicals import *
from chemicals.identifiers import pubchem_db

try:
    Chemical('asdfasdfsadfasd')
except:
    pass


def get_radius_gyration(lines, MW, linear=False):
    '''Given a molecular weight and input lines from psi4
    calculate the radius of gyration in units of m.
    '''
    rotational_constants_line = 'Rotational constants: '
    for l in lines:
        if rotational_constants_line in l:
            l = l.strip()
            segments = l.replace(rotational_constants_line, '').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').split(' ')
            #print(l)
            #print(segments)
            assert segments[-1] == '[cm^-1]'
            _, _, A, _, _, B, _, _, C, _ = segments
            try:
                A, B, C = float(A), float(B), float(C)
                A, B, C = A*c*100.0, B*c*100.0, C*c*100.0
                A, B, C = (h/(8*pi**2)/i for i in (A, B, C))
            except:
                if '*' in A:
                    A, B = float(B), float(C)
                    A, B = A*c*100.0, B*c*100.0
                    A, B = (h/(8*pi**2)/i for i in (A, B))
                    linear = True
                else:
                    raise NotImplementedError('DEBUG')

            if not linear:
                return radius_of_gyration(A=A, B=B, C=C, MW=MW)
            else:
                return radius_of_gyration(A=A, B=B, C=0, MW=MW, planar=True)
    raise ValueError("Rotational constants not found")

def get_dipole_moment(lines):
    expect_line = '  Dipole Moment: [D]'
    for i, l in enumerate(lines):
#         print(l)
        if l.startswith(expect_line):
            l = lines[i+1]
            dipole_moment = l.split(' ')[-1]
            dipole_moment = float(dipole_moment)
#             print(dipole_moment)
            return dipole_moment
    raise ValueError("Dipole moment not found")


def get_is_linear(lines):
    # Never seems to be linear in the output
    expect_line = '    Rotational Cp'
    for i, l in enumerate(lines):
        if l.startswith(expect_line):
            return '12.472' not in l
    raise ValueError("Cp not found")

def get_total_Cp(lines):
    expect_line = '  Total Cp'
    for i, l in enumerate(lines):
        if l.startswith(expect_line):
            l = l.replace(expect_line, '')
            while '  ' in l:
                l = l.replace('  ', ' ')
            segments = l.split(' ')
            Cp = segments[4]
            Cp = float(Cp)
            return Cp
    raise ValueError("Total Cp not found")

def get_vibration_freqs(lines):
    expect_line = '  Freq [cm^-1]'
    freqs = []
    for i, l in enumerate(lines):
        if l.startswith(expect_line):
            #print(l)
            l = l.strip().replace(expect_line, '').replace(expect_line, '').replace('Freq', '').replace('[cm^-1]', '')
            while '  ' in l:
                l = l.replace('  ', ' ')
            segments = l.split(' ')
            #print(segments)
            for s in segments:
                if s:
                    if s.endswith('i'):
                        # Probably a bug and the results should be discarded
                        continue
                    freqs.append(float(s))
#             print(segments)
    if not freqs:
        raise ValueError("Frequencies not found")
    return freqs


working_folder = 'psi4'
to_process = []
for f in os.listdir(working_folder):
    if f.endswith('.out'):
        to_process.append(os.path.join(working_folder, f))

dipoles = {}
radius_gyrations = {}
linears = {}
Cps = {}
frequencies = {}

def get_data(f):
    lines = open(f).readlines()
    CAS = f.replace(working_folder, '').replace('.out', '').replace('/', '')

    failure_OK = False
    for l in lines[-5:]:
        if 'Psi4 encountered an error' in l:
            failure_OK = True

    try:
        dipole = get_dipole_moment(lines)
        dipoles[CAS] = dipole
    except Exception as e:
        if not failure_OK:
            raise e

    try:
        Cp = get_total_Cp(lines)
        Cps[CAS] = Cp
    except Exception as e:
        if not failure_OK:
            print(CAS)
            raise e

    try:
        linear = get_is_linear(lines)
        linears[CAS] = linear
    except Exception as e:
        if not failure_OK:
            print(CAS)
            raise e



    try:
        freqs = get_vibration_freqs(lines)
        frequencies[CAS] = freqs
    except ValueError as e:
        if not failure_OK:
            print(CAS)
            raise e

    MW = pubchem_db.search_CAS(CAS).MW

    try:
        rgyr = get_radius_gyration(lines, MW, linear=linear)
        radius_gyrations[CAS] = rgyr
    except Exception as e:
        if not failure_OK:
            print(CAS)
            raise e

for f in to_process:
    get_data(f)

#print(dipoles)
print(f'Found {len(dipoles)} dipoles' )
print(f'Found {len(linears)} linears' )
print(f'Found {len(frequencies)} frequencies' )
print(f'Found {len(radius_gyrations)} radius_gyrations' )

#for k, v in linears.items():
    #if v:
        #print(k)

#print(frequencies)


def write_dipole_file():
    dump_oder = ['dipole_moment',]
    keys = ['CAS', 'dipole_moment']
    lines = ['\t'.join(keys) + '\n']
    for CAS, dipole in dipoles.items():
        values = [dipole]
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
        lines.append(to_write)

    f = open('../../chemicals/Misc/psi4_dipoles.tsv', 'w')
    f.writelines(lines)
    f.close()
    return True

def write_rg_file():
    dump_oder = ['RG',]
    keys = ['CAS', 'RG']
    lines = ['\t'.join(keys) + '\n']
    for CAS, rg in radius_gyrations.items():
        values = [rg]
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
        lines.append(to_write)

    f = open('../../chemicals/Misc/psi4_radius_of_gyrations.tsv', 'w')
    f.writelines(lines)
    f.close()
    return True

def write_linear_file():
    dump_oder = ['linear',]
    keys = ['CAS', 'linear']
    lines = ['\t'.join(keys) + '\n']
    for CAS, linear in linears.items():
        values = [linear]
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
        lines.append(to_write)

    f = open('../../chemicals/Misc/psi4_linear.tsv', 'w')
    f.writelines(lines)
    f.close()
    return True

adjusted_Ts = {}
unadjusted_Ts = {}
def write_frequencies():
    for CAS, cminvs in frequencies.items():
        adjusted_Ts[CAS] = [vibration_frequency_cm_to_characteristic_temperature(f, scale=0.9365) for f in cminvs]
        unadjusted_Ts[CAS] = [vibration_frequency_cm_to_characteristic_temperature(f, scale=1) for f in cminvs]

    separators = (',', ':')
    json.dump(adjusted_Ts, open('../../chemicals/Heat Capacity/psi4_adjusted_characteristic_temperatures.json', 'w'),
              sort_keys=True, indent=2, separators=separators)
    json.dump(unadjusted_Ts, open('../../chemicals/Heat Capacity/psi4_unadjusted_characteristic_temperatures.json', 'w'),
              sort_keys=True, indent=2, separators=separators)

write_dipole_file()
write_rg_file()
write_frequencies()
write_linear_file()
#print(Cps)
#def print_Cp_diffs():
    #for CAS, temps in unadjusted_Ts.items():
        #Cp_calc = Cpg_statistical_mechanics(298.15, temps, linear=linears[CAS])
        #print(Cp_calc/Cps[CAS], Cp_calc, Cps[CAS], linears[CAS], temps, CAS)

#print_Cp_diffs()
