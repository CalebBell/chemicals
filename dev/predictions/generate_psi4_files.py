# -*- coding: utf-8 -*-
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
SOFTWARE.'''

import sys, os, argparse
from chemicals import *
from thermo import *
from rdkit import Chem
from joblib import Parallel, delayed
from subprocess import Popen, PIPE


from fluids.constants import c, h, pi


def get_radius_gyration(lines, MW, linear=False):
    '''Given a molecular weight and input lines from psi4
    calculate the radius of gyration in units of m.
    '''
    rotational_constants_line = '  Rotational constants: '
    for l in lines:
        if l.startswith(rotational_constants_line):
            segments = l.replace(rotational_constants_line, '').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').split(' ')
            assert segments[-1] == '[cm^-1]'
            _, _, A, _, _, B, _, _, C, _ = segments
            A, B, C = float(A), float(B), float(C)
            A, B, C = A*c*100.0, B*c*100.0, C*c*100.0
            A, B, C = [h/(8*pi**2)/i for i in (A, B, C)]
            
            if not linear:
                return radius_of_gyration(A=A, B=B, C=C, MW=MW)
            else:
                return radius_of_gyration(A=A, B=B, MW=MW, planar=True)
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
            l = l.replace(expect_line, '')
            while '  ' in l:
                l = l.replace('  ', ' ')
            segments = l.split(' ')
            for s in segments:
                if s:
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
    lines = open(f, 'r').readlines()
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
        linear = get_is_linear(lines)
        linears[CAS] = linear
    except Exception as e:
        if not failure_OK:
            print(CAS)
            raise e

for f in to_process:
    get_data(f)
print(dipoles)
print(f'Found {len(dipoles)} dipoles' )

