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

import os
from random import shuffle
from subprocess import PIPE, Popen

from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import AllChem
from thermo import *

from chemicals import *
from chemicals.identifiers import pubchem_db


def is_psi4_fast(metadata):
    try:
        atoms = nested_formula_parser(metadata.formula)
    except:
        # Failing on isotopes
        return False
    atom_count = sum(v for e, v in atoms.items() if e != 'H')
    is_small = atom_count < 4
    #is_small = atom_count == 3
    return is_small

def rdkit_atom_position_as_xyz(rdkit_atom):
    # TODO use better formatting
    line = "{} {} {} {}"
    conformer = rdkit_atom.GetOwningMol().GetConformer()
    position = conformer.GetAtomPosition(rdkit_atom.GetIdx())

    return f"{rdkit_atom.GetSymbol()} {position.x} {position.y} {position.z}\n\n"
    return line

def mol_with_Hs(mol):
    # Generaate the molecule with Hs
    if mol.GetNumConformers() <= 1:
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        except:
            # Oh well
            pass
    else:
        mol = Chem.AddHs(mol)
    return mol

# force loading all data
try:
    Chemical('asdfasdfsadfasd')
except:
    pass

def write_input(mol, dest):
    mol_Hs = mol_with_Hs(mol)

    molxyz = """
#! psi4input\n\n

molecule themol
{\n0 1

"""


    remainder = """
set T 298.15
set P 101325
set basis 6-31G**
#set basis 6-311++G(3df,3pd)
set reference uhf
optimize( "mp2" )
#energy( "mp2" )
freq('mp2',molecule=themol)
    """
    atoms = mol_Hs.GetAtoms()

    for atom in atoms:
        molxyz += rdkit_atom_position_as_xyz(atom)
    molxyz += "}\n"

    f = open(dest, 'w')
    f.write(molxyz)
    f.write(remainder)
    f.close()
    return True


def call_psi4(CAS):
        file_in = os.path.join(working_folder, "%s.dat" %(CAS))
        file_out = os.path.join(working_folder, "%s.out" %(CAS))

        run = True
        if os.path.exists(file_out):
            output_file = open(file_out)
            output_lines = output_file.readlines()
            output_file.close()
            for l in output_lines[-5:]:
                if '*** Psi4 exiting successfully.' in l:
                    run = False
                elif 'Psi4 encountered an error' in l:
                    run = False
                if not run:
                    break
            if run:
                os.remove(file_out)
        if run:
            process = Popen(['psi4',file_in, file_out, '-n', '1'], stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()

        if os.path.exists(file_out):
            output_file = open(file_out)
            output_lines = output_file.readlines()
            output_file.close()

            # Remove nonreproducible lines
            output_lines = [i for i in output_lines if 'Psi4 started on' not in i]
            output_lines = [i for i in output_lines if 'Process ID' not in i]
            output_lines = [i for i in output_lines if 'Psi4 stopped on' not in i]
            output_lines = [i for i in output_lines if 'Psi4 wall time' not in i]
            output_lines = [i for i in output_lines if '*** at' not in i]
            output_lines = [i for i in output_lines if '*** tstop()' not in i]
            output_lines = [i for i in output_lines if '	user time' not in i]
            output_lines = [i for i in output_lines if '	system time' not in i]
            output_lines = [i for i in output_lines if '	total time' not in i]


            output_file = open(file_out, 'w')
            output_file.writelines(output_lines)
            output_file.close()
        return True


working_folder = 'psi4'

to_process = []
for metadata in pubchem_db.CAS_index.values():
    if is_psi4_fast(metadata):
        CAS = metadata.CASs
        c = Chemical(CAS)
        mol = c.rdkitmol
        #print(c, c.atoms)
        if mol is None:
            continue
        dest = os.path.join(working_folder, "%s.dat" %(metadata.CASs))
        write_input(mol, dest)
        to_process.append(metadata.CASs)
shuffle(to_process)

        #psi4 50-00-0.dat 50-00-0.out -n 1
retVals = Parallel(n_jobs=7, verbose=100, backend='threading')(delayed(call_psi4)(CAS) for CAS in sorted(to_process))
