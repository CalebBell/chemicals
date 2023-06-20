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
import xml.etree.ElementTree as ET

folder = os.path.join(os.path.dirname(__file__), '..','..', 'chemicals', 'Misc')

tree = ET.parse(os.path.join(folder, 'ChemSep8.30.xml'))
root = tree.getroot()

radius_of_gyrations = {}
for child in root:
    CAS, name, smiles, formula, rgyr = None, None, None, None, None
    for i in child:
        tag = i.tag
        if CAS is None and tag == 'CAS':
            CAS = i.attrib['value']
        if rgyr is None and tag == 'RadiusOfGyration':
            rgyr = i.attrib['value']
        if rgyr is not None and CAS is not None:
            radius_of_gyrations[CAS] = float(rgyr)

# sort it so if the chemsep file gets updated the file has minimal changes
radius_of_gyrations = {k: radius_of_gyrations[k] for k in sorted(radius_of_gyrations.keys())}

dump_oder = ['RG',]
keys = ['CAS', 'RG']
lines = ['\t'.join(keys) + '\n']
for CAS, RG in radius_of_gyrations.items():
    values = [RG]
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

f = open('../../chemicals/Misc/chemsep_radius_of_gyrations.tsv', 'w')
f.writelines(lines)
f.close()
