# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell
<Caleb.Andrew.Bell@gmail.com> 
Copyright (C) 2020 Yoel Rene Cortes-Pena
<yoelcortes@gmail.com>

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

This module contains a complete periodic table, routines for working with
chemical formulas, computing molecular weight, computing mass fractions and 
atom fractions, and assorted other tasks.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Periodic Table and Elements
---------------------------
.. autodata:: chemicals.elements.periodic_table
.. autoclass:: chemicals.elements.Element
.. autoclass:: chemicals.elements.PeriodicTable
          
Working with Formulas
---------------------
.. autofunction:: chemicals.elements.simple_formula_parser
.. autofunction:: chemicals.elements.nested_formula_parser
.. autofunction:: chemicals.elements.charge_from_formula
.. autofunction:: chemicals.elements.serialize_formula
.. autofunction:: chemicals.elements.atoms_to_Hill

Working with Parsed Formulas
----------------------------
.. autofunction:: chemicals.elements.molecular_weight
.. autofunction:: chemicals.elements.similarity_variable
.. autofunction:: chemicals.elements.atom_fractions
.. autofunction:: chemicals.elements.mass_fractions
.. autofunction:: chemicals.elements.mixture_atomic_composition
.. autofunction:: chemicals.elements.mixture_atomic_composition_ordered
.. autofunction:: chemicals.elements.atom_matrix

"""

__all__ = ['PeriodicTable', 'molecular_weight', 'mass_fractions', 
           'atom_fractions','mixture_atomic_composition', 'atom_matrix',
           'similarity_variable', 'atoms_to_Hill', 
           'simple_formula_parser', 'nested_formula_parser', 'CAS_by_number',
           'periods', 'groups',  'homonuclear_elements',
           'blocks', 'homonuclear_elemental_gases', 'charge_from_formula',
           'serialize_formula', 'mixture_atomic_composition_ordered',
           'periodic_table']
import re

CAS_by_number_standard = ['1333-74-0', '7440-59-7', '7439-93-2', '7440-41-7', '7440-42-8', '7440-44-0', '7727-37-9', '7782-44-7', '7782-41-4', '7440-01-9', '7440-23-5', '7439-95-4', '7429-90-5', '7440-21-3', '7723-14-0', '7704-34-9', '7782-50-5', '7440-37-1', '7440-09-7', '7440-70-2', '7440-20-2', '7440-32-6', '7440-62-2', '7440-47-3', '7439-96-5', '7439-89-6', '7440-48-4', '7440-02-0', '7440-50-8', '7440-66-6', '7440-55-3', '7440-56-4', '7440-38-2', '7782-49-2', '10097-32-2', '7439-90-9', '7440-17-7', '7440-24-6', '7440-65-5', '7440-67-7', '7440-03-1', '7439-98-7', '7440-26-8', '7440-18-8', '7440-16-6', '7440-05-3', '7440-22-4', '7440-43-9', '7440-74-6', '7440-31-5', '7440-36-0', '13494-80-9', '7553-56-2', '7440-63-3', '7440-46-2', '7440-39-3', '7439-91-0', '7440-45-1', '7440-10-0', '7440-00-8', '7440-12-2', '7440-19-9', '7440-53-1', '7440-54-2', '7440-27-9', '7429-91-6', '7440-60-0', '7440-52-0', '7440-30-4', '7440-64-4', '7439-94-3', '7440-58-6', '7440-25-7', '7440-33-7', '7440-15-5', '7440-04-2', '7439-88-5', '7440-06-4', '7440-57-5', '7439-97-6', '7440-28-0', '7439-92-1', '7440-69-9', '7440-08-6', '7440-68-8', '10043-92-2', '7440-73-5', '7440-14-4', '7440-34-8', '7440-29-1', '7440-13-3', '7440-61-1', '7439-99-8', '7440-07-5', '7440-35-9', '7440-51-9', '7440-40-6', '7440-71-3', '7429-92-7', '7440-72-4', '7440-11-1', '10028-14-5', '22537-19-5', '53850-36-5', '53850-35-4', '54038-81-2', '54037-14-8', '54037-57-9', '54038-01-6', '54083-77-1', '54386-24-2', '54084-26-3', '54084-70-7', '54085-16-4', '54085-64-2', '54100-71-9', '54101-14-3', '54144-19-3']
CAS_by_number = CAS_by_number_standard.copy()
'''CAS numbers of the elements, indexed by atomic numbers off-by-one up to 118.'''

periods = [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
'''Periods of the elements, indexed by atomic numbers off-by-one up to 118.'''

groups = [1, 18, 1, 2, 13, 14, 15, 16, 17, 18, 1, 2, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
'''Groups of the elements, indexed by atomic numbers off-by-one up to 118.
Lanthanides and Actinides are set to None.'''

s_block = [1, 2, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88]
d_block = [*range(21, 31), *range(39, 49), *range(72, 81), *range(104, 113)]
f_block = [*range(57, 72), *range(89, 104)] # 57, 89 are sometimes placed in the d block
p_block = [*range(5, 11), *range(13, 19), *range(31, 37), *range(49, 55), *range(81, 87), *range(113, 119)]
blocks = {'s': s_block, 'd': d_block, 'f': f_block, 'p': p_block}
'''Blocks of the elements, stored in a dictionary with four keys and lists.
Indexed by atomic numbers off-by-one up to 118.'''

InChI_keys = ['YZCKVEUIGOORGS-UHFFFAOYSA-N', 'SWQJXJOGLNCZEY-UHFFFAOYSA-N', 'WHXSMMKQMYFTQS-UHFFFAOYSA-N', 'ATBAMAFKBVZNFJ-UHFFFAOYSA-N', 'ZOXJGFHDIHLPTG-UHFFFAOYSA-N', 'OKTJSMMVPCPJKN-UHFFFAOYSA-N', 'QJGQUHMNIGDVPM-UHFFFAOYSA-N', 'QVGXLLKOCUKJST-UHFFFAOYSA-N', 'YCKRFDGAMUMZLT-UHFFFAOYSA-N', 'GKAOGPIIYCISHV-UHFFFAOYSA-N', 'KEAYESYHFKHZAL-UHFFFAOYSA-N', 'FYYHWMGAXLPEAU-UHFFFAOYSA-N', 'XAGFODPZIPBFFR-UHFFFAOYSA-N', 'XUIMIQQOPSSXEZ-UHFFFAOYSA-N', 'OAICVXFJPJFONN-UHFFFAOYSA-N', 'NINIDFKCEFEMDL-UHFFFAOYSA-N', 'ZAMOUSCENKQFHK-UHFFFAOYSA-N', 'XKRFYHLGVUSROY-UHFFFAOYSA-N', 'ZLMJMSJWJFRBEC-UHFFFAOYSA-N', 'OYPRJOBELJOOCE-UHFFFAOYSA-N', 'SIXSYDAISGFNSX-UHFFFAOYSA-N', 'RTAQQCXQSZGOHL-UHFFFAOYSA-N', 'LEONUFNNVUYDNQ-UHFFFAOYSA-N', 'VYZAMTAEIAYCRO-UHFFFAOYSA-N', 'PWHULOQIROXLJO-UHFFFAOYSA-N', 'XEEYBQQBJWHFJM-UHFFFAOYSA-N', 'GUTLYIVDDKVIGB-UHFFFAOYSA-N', 'PXHVJJICTQNCMI-UHFFFAOYSA-N', 'RYGMFSIKBFXOCR-UHFFFAOYSA-N', 'HCHKCACWOHOZIP-UHFFFAOYSA-N', 'GYHNNYVSQQEPJS-UHFFFAOYSA-N', 'GNPVGFCGXDBREM-UHFFFAOYSA-N', 'RQNWIZPPADIBDY-UHFFFAOYSA-N', 'BUGBHKTXTAQXES-UHFFFAOYSA-N', 'WKBOTKDWSSQWDR-UHFFFAOYSA-N', 'DNNSSWSSYDEUBZ-UHFFFAOYSA-N', 'IGLNJRXAVVLDKE-UHFFFAOYSA-N', 'CIOAGBVUUVVLOB-UHFFFAOYSA-N', 'VWQVUPCCIRVNHF-UHFFFAOYSA-N', 'QCWXUUIWCKQGHC-UHFFFAOYSA-N', 'GUCVJGMIXFAOAE-UHFFFAOYSA-N', 'ZOKXTWBITQBERF-UHFFFAOYSA-N', 'GKLVYJBZJHMRIY-UHFFFAOYSA-N', 'KJTLSVCANCCWHF-UHFFFAOYSA-N', 'MHOVAHRLVXNVSD-UHFFFAOYSA-N', 'KDLHZDBZIXYQEI-UHFFFAOYSA-N', 'BQCADISMDOOEFD-UHFFFAOYSA-N', 'BDOSMKKIYDKNTQ-UHFFFAOYSA-N', 'APFVFJFRJDLVQX-UHFFFAOYSA-N', 'ATJFFYVFTNAWJD-UHFFFAOYSA-N', 'WATWJIUSRGPENY-UHFFFAOYSA-N', 'PORWMNRCUJJQNO-UHFFFAOYSA-N', 'ZCYVEMRRCGMTRW-UHFFFAOYSA-N', 'FHNFHKCVQCLJFQ-UHFFFAOYSA-N', 'TVFDJXOCXUVLDH-UHFFFAOYSA-N', 'DSAJWYNOEDNPEQ-UHFFFAOYSA-N', 'FZLIPJUXYLNCLC-UHFFFAOYSA-N', 'GWXLDORMOJMVQZ-UHFFFAOYSA-N', 'PUDIUYLPXJFUGB-UHFFFAOYSA-N', 'QEFYFXOXNSNQGX-UHFFFAOYSA-N', 'VQMWBBYLQSCNPO-UHFFFAOYSA-N', 'KZUNJOHGWZRPMI-UHFFFAOYSA-N', 'OGPBJKLSAFTDLK-UHFFFAOYSA-N', 'UIWYJDYFSGRHKR-UHFFFAOYSA-N', 'GZCRRIHWUXGPOV-UHFFFAOYSA-N', 'KBQHZAAAGSGFKK-UHFFFAOYSA-N', 'KJZYNXUDTRRSPN-UHFFFAOYSA-N', 'UYAHIZSMUZPPFV-UHFFFAOYSA-N', 'FRNOGLGSGLTDKL-UHFFFAOYSA-N', 'NAWDYIZEMPQZHO-UHFFFAOYSA-N', 'OHSVLFRHMCKCQY-UHFFFAOYSA-N', 'VBJZVLUMGGDVMO-UHFFFAOYSA-N', 'GUVRBAGPIYLISA-UHFFFAOYSA-N', 'WFKWXMTUELFFGS-UHFFFAOYSA-N', 'WUAPFZMCVAUBPE-UHFFFAOYSA-N', 'SYQBFIAQOQZEGI-UHFFFAOYSA-N', 'GKOZUEZYRPOHIO-UHFFFAOYSA-N', 'BASFCYQUMIYNBI-UHFFFAOYSA-N', 'PCHJSUWPFVWCPO-UHFFFAOYSA-N', 'QSHDDOUJBYECFT-UHFFFAOYSA-N', 'BKVIYDNLLOSFOA-UHFFFAOYSA-N', 'WABPQHHGFIMREM-UHFFFAOYSA-N', 'JCXGWMGPZLAOME-UHFFFAOYSA-N', 'HZEBHPIOVYHPMT-UHFFFAOYSA-N', 'RYXHOMYVWAEKHL-UHFFFAOYSA-N', 'SYUHGPGVQRZVTB-UHFFFAOYSA-N', 'KLMCZVJOEAUDNE-UHFFFAOYSA-N', 'HCWPIIXVSYCSAN-UHFFFAOYSA-N', 'QQINRWTZWGJFDB-UHFFFAOYSA-N', 'ZSLUVFAKFWKJRC-UHFFFAOYSA-N', 'XLROVYAPLOFLNU-UHFFFAOYSA-N', 'JFALSRSLKYAFGM-UHFFFAOYSA-N', 'LFNLGNPSGWYGGD-UHFFFAOYSA-N', 'OYEHPCDNVJXUIW-UHFFFAOYSA-N', 'LXQXZNRPTYVCNG-UHFFFAOYSA-N', 'NIWWFAAXEMMFMS-UHFFFAOYSA-N', 'PWVKJRSRVJTHTR-UHFFFAOYSA-N', 'HGLDOAKPQXAFKI-UHFFFAOYSA-N', 'CKBRQZNRCSJHFT-UHFFFAOYSA-N', 'MIORUQGGZCBUGO-UHFFFAOYSA-N', 'MQVSLOYRCXQRPM-UHFFFAOYSA-N', 'ORQBXQOJMQIAOY-UHFFFAOYSA-N', 'CNQCVBJFEGMYDW-UHFFFAOYSA-N', 'YGPLJIIQQIDVFJ-UHFFFAOYSA-N', 'PUKKTGLVJQVIOF-UHFFFAOYSA-N', 'VAOUCABZIBBBJH-UHFFFAOYSA-N', 'INOXRQQPOOCQPH-UHFFFAOYSA-N', 'OBDWMWVOVYJOMI-UHFFFAOYSA-N', 'VAJSJTKWMRUWBF-UHFFFAOYSA-N', 'NCBMSFCPDGXTHD-UHFFFAOYSA-N', 'LJROPTGWFUZRDB-UHFFFAOYSA-N', 'NOTIIDSZELDPOP-UHFFFAOYSA-N', 'KUGNSLWRKGRKGS-UHFFFAOYSA-N', 'WIHJCBVMYKIGOT-UHFFFAOYSA-N', 'QDXZEHQJHSHEQF-UHFFFAOYSA-N', 'ONFASNXETZOODS-UHFFFAOYSA-N', 'INMSAURDCVBGHH-UHFFFAOYSA-N', 'GOANEQIZDYDFCO-UHFFFAOYSA-N']
# Big problem: Atoms like N2, O2 point to only the singlet

homonuclear_elements = {1, 7, 8, 9, 17, 35, 53} # includes Br2, I2
# homonuclear_elemental_gases = [periodic_table[i] for i in ['hydrogen', 'nitrogen', 'oxygen', 'fluorine', 'chlorine', 'bromine', 'iodine']]
homonuclear_elemental_gases = [1, 7, 8, 9, 17] # 35, 53
homonuclear_elemental_singlets_CASs = ["12385-13-6", "17778-88-0", "17778-80-2", "14762-94-8", "22537-15-1"]
for i, CAS in zip(homonuclear_elemental_gases, homonuclear_elemental_singlets_CASs):
    CAS_by_number[i-1] = CAS

cids = [783, 23987, 3028194, 5460467, 5462311, 5462310, 222, 962, 14917, 23935, 5360545, 5462224, 5359268, 5461123, 24404, 402, 313, 23968, 5462222, 5460341, 23952, 23963, 23990, 23976, 23930, 23925, 104730, 935, 23978, 23994, 5360835, 6326954, 5359596, 6326970, 260, 5416, 5357696, 5359327, 23993, 23995, 23936, 23932, 23957, 23950, 23948, 23938, 23954, 23973, 5359967, 5352426, 5354495, 6327182, 24841, 23991, 5354618, 5355457, 23926, 23974, 23942, 23934, 23944, 23951, 23981, 23982, 23958, 23912, 23988, 23980, 23961, 23992, 23929, 23986, 23956, 23964, 23947, 23937, 23924, 23939, 23985, 23931, 5359464, 5352425, 5359367, 6328143, 5460479, 24857, 6328145, 6328144, 23965, 23960, 23945, 23989, 23933, 23940, 23966, 23979, 23971, 23997, 23913, 23998, 23943, 24822, 31192, 56951715, 56951718, 56951717, 56951713, 56951714, 56951716, None, None, None, None, None, None, None, None, None]

'''
liquids = ['Hg', 'Br']
gases = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'F', 'Cl', 'H', 'O', 'N']
phases = ['s' for i in periodic_table]
for i, ele in enumerate(periodic_table):
    if ele.symbol in liquids:
        phases[i] = 'l'
    elif ele.symbol in gases:
        phases[i] = 'g'
'''
phases = ['g', 'g', 's', 's', 's', 's', 'g', 'g', 'g', 'g', 's', 's', 's', 's', 's', 's', 'g', 'g', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 'l', 'g', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 'g', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 'l', 's', 's', 's', 's', 's', 'g', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's', 's']

'''# From CRC Table:
from chemicals.heat_capacity import CRC_standard_data
S0s = []
Hfs = []

other_values = {
    '7440-31-5': (0, 51.18), # Tin (white) instead of gray from CODATA
    '10097-32-2': (0, 152.2) # bromine as a liquid
}
for ele in periodic_table:
    try:
        CAS = ele.CAS_standard
        if CAS == '7440-44-0': # carbon -> graphite (which is the standard state)
            CAS = '7782-42-5'
        S0c = CRC_standard_data.at[CAS, 'S0s']
        S0l = CRC_standard_data.at[CAS, 'S0l']
        S0g = CRC_standard_data.at[CAS, 'S0g']
        Hfs = CRC_standard_data.at[CAS, 'Hfs']
        Hfl = CRC_standard_data.at[CAS, 'Hfl']
        Hfg = CRC_standard_data.at[CAS, 'Hfg']        
        if ele.phase == 's':
            S0 = S0c
            Hf = Hfs
        elif ele.phase == 'l':
            S0 = S0l
            Hf = Hfl
        elif ele.phase == 'g':
            S0 = S0g
            Hf = Hfg
        if isnan(S0):
            S0 = None
        if isnan(Hf):
            Hf = None
        if CAS in other_values:
            Hf, S0 = other_values[CAS]
            
        S0s.append(S0)
        Hfs.append(Hf)

    except Exception as e:
        S0s.append(None)
        Hfs.append(None)

for Hf, ele in zip(Hfs, periodic_table):
    if Hf != 0 and Hf is not None:
        print(Hf, ele.name, ele.CAS)
'''

# Note that atoms like Br2, I2, O2, N2 have values of S for two atoms not one
S0s = [130.7, 126.2, 29.1, 9.5, 5.9, 5.7, 191.6, 205.2, 202.8, 146.3, 51.3, 32.7, 28.3, 18.8, 41.1, 32.1, 223.1, 154.8, 64.7, 41.6, 34.6, 30.7, 28.9, 23.8, 32.0, 27.3, 30.0, 29.9, 33.2, 41.6, 40.8, 31.1, 35.1, 42.4, 152.2, 164.1, 76.8, 55.0, 44.4, 39.0, 36.4, 28.7, None, 28.5, 31.5, 37.6, 42.6, 51.8, 57.8, 51.18, 45.7, 49.7, 116.1, 169.7, 85.2, 62.5, 56.9, 72.0, 73.2, 71.5, None, 69.6, 77.8, 68.1, 73.2, 75.6, 75.3, 73.2, 74.0, 59.9, 51.0, 43.6, 41.5, 32.6, 36.9, 32.6, 35.5, 41.6, 47.4, 75.9, 64.2, 64.8, 56.7, None, None, 176.2, 95.4, 71.0, 56.5, 51.8, 51.9, 50.2, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
Hfs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]

class PeriodicTable:
    """Periodic Table object for use in dealing with elements.

    As there is only one periodic table of elements, this is automatically
    initialized into the object `periodic_table`; there is no need to 
    construct a new instance of this class.

    Parameters
    ----------
    elements : list[Element]
        List of Element objects

    Notes
    -----
    Can be checked to sese if an element in in this, can be iterated over,
    and as a current length of 118 elements.

    See Also
    --------
    periodic_table
    Element

    References
    ----------
    .. [1] N M O'Boyle, M Banck, C A James, C Morley, T Vandermeersch, and
       G R Hutchison. "Open Babel: An open chemical toolbox." J. Cheminf.
       (2011), 3, 33. DOI:10.1186/1758-2946-3-33
    """
    __slots__ = ('_number_to_elements', '_symbol_to_elements',
                 '_name_to_elements', '_CAS_to_elements', '_indexes')
    def __init__(self, elements):
        #: Dictionary lookup of number(int) -> Element;
        #: also has number(str) -> Element for convenience.
        self._number_to_elements = number_to_elements = {}        
        #: Dictionary lookup of symbol(str) -> Element.
        self._symbol_to_elements = symbol_to_elements = {}        
        #: Dictionary lookup of name(str) -> Element;
        #: also has name(str.lower()) -> Element for convenience.
        self._name_to_elements = name_to_elements = {}
        #: Dictionary lookup of CAS(str) -> Element.
        self._CAS_to_elements = CAS_to_elements = {}
    
        for ele in elements:
            number_to_elements[ele.number] = ele
            number_to_elements[str(ele.number)] = ele
            symbol_to_elements[ele.symbol] = ele
            name_to_elements[ele.name] = ele
            name_to_elements[ele.name.lower()] = ele
            CAS_to_elements[ele.CAS] = ele
            
        self._indexes = (symbol_to_elements, number_to_elements, 
                         name_to_elements, CAS_to_elements)

    def __contains__(self, key):
        for i in self._indexes:
            if key in i: return True
        return False
    
    def __len__(self):
        return 118

    def __iter__(self):
        return iter([self._number_to_elements[i] for i in range(1,119)])

    def __getitem__(self, key):
        for i in self._indexes:
            if key in i: return i[key]
        raise KeyError("'%s' is not in the periodic table" %(key))

    def __getattr__(self, key):
        for i in self._indexes:
            if key in i: return i[key]
        raise AttributeError("'%s' is not in the periodic table" %(key))


class Element:
    """Class for storing data on chemical elements. Supports most common
    properties. If a property is not available, it is set to None.
    
    The elements are created automatically and should be accessed via the 
    `periodic_table` interface.

    Attributes
    ----------
    number : int
        Atomic number
    name : str
        name
    symbol : str
        Elemental symbol
    MW : float
        Molecular weight
    CAS : str
        CAS number
    period : str
        Period in the periodic table
    group : str
        Group in the periodic table
    block : str
        Block in the periodic table
    AReneg : float
        Allred and Rochow electronegativity
    rcov : float
        Covalent radius, [Angstrom]
    rvdw : float
        Van der Waals radius, [Angstrom]
    maxbonds : float
        Maximum valence of a bond with this element
    elneg : float
        Pauling electronegativity
    ionization : float
        Ionization potential, [eV]
    ionization : float
        elaffinity affinity, [eV]
    protons : int
        Number of protons
    electrons : int
        Number of electrons of the element in the ground state
    InChI : str
        Standard InChI string of the element
    InChI_key : str
        25-character hash of the compound's InChI.
    smiles : str
        Standard smiles string of the element
    PubChem : int
        PubChem Compound identifier (CID) of the chemical
    phase : str
        Standard state at 1 atm and 298.15 K
    Hf : float
        Enthalpy of formation of the element in its standard state (0 by
        definition), [J/mol]
    S0 : float
        Standard absolute entropy of the element in its standard state (1 bar,
        298.15 K), [J/mol/K]
    """
    __slots__ = ['number', 'symbol', 'name', 'CAS', 'MW', 'AReneg', 'rcov',
                 'rvdw', 'maxbonds', 'elneg', 'ionization', 'elaffinity',
                 'period', 'group', 
                 'InChI_key', 'PubChem', 'phase', 'Hf', 'S0']
    
    def __repr__(self):
        return "<Element %s (%s), number %d, MW=%s>" %(self.name, self.symbol, self.number, self.MW)
    
    def __init__(self, number, symbol, name, MW, CAS, AReneg, rcov, rvdw,
                 maxbonds, elneg, ionization, elaffinity, period, group,
                 PubChem, phase, Hf, S0, InChI_key=None):
        self.number = number
        self.symbol = symbol
        self.name = name
        self.MW = MW
        self.CAS = CAS

        self.period = period
        self.group = group

        self.AReneg = AReneg
        self.rcov = rcov
        self.rvdw = rvdw
        self.maxbonds = maxbonds
        self.elneg = elneg
        self.ionization = ionization
        self.elaffinity = elaffinity
        
        self.InChI_key = InChI_key
        self.PubChem = PubChem
        
        self.phase = phase
        self.S0 = S0
        self.Hf = Hf

    @property
    def CAS_standard(self):
        r'''CAS number of the compound of the element used as a standard state
        ; i.e. the typically diatomic molecules hydrogen, nitrogen, oxygen, 
        fluorine, and chlorine, have different CAS numbers for the monoatomic
        form and the diatomic form. This method returns the conventionally used
        CAS number.
        '''
        return CAS_by_number_standard[self.number-1]
    
    @property
    def protons(self): 
        r'''The number of protons of the element.'''
        return self.number
    
    @property
    def electrons(self): 
        r'''The number of electrons of the element.'''
        return self.number

    @property
    def smiles(self): 
        r'''The SMILES identification string of the element.'''
        return '[' + self.symbol + ']'

    @property
    def InChI(self): 
        r'''The InChI identifier of the element. One of 's', 'd', 'f', or 'p'.'''
        return self.symbol # 'InChI=1S/' +
    
    @property
    def block(self):
        r'''Which block of the periodic table the element is in.'''
        for k, v in blocks.items():
            if self.number in v:
                return k


# New file format; same data, coverted to Python lists for convenience with a few regular expressions
# https://github.com/openbabel/openbabel/blob/master/src/elementtable.h
openbabel_element_data = [
[  1, "H", 2.20, 0.31, 0.31, 1.10,  1,    1.00794, 2.20, 13.5984, 0.75420375, 0.75, 0.75, 0.75, "Hydrogen"],
[  2, "He", 0.00, 0.28, 0.28, 1.40,  0,   4.002602, 0.00, 24.5874,          0, 0.85, 1.00, 1.00, "Helium"],
[  3, "Li", 0.97, 1.28, 1.28, 1.81,  1,      6.941, 0.98,  5.3917,   0.618049, 0.80, 0.50, 1.00, "Lithium"],
[  4, "Be", 1.47, 0.96, 0.96, 1.53,  2,   9.012182, 1.57,  9.3227,          0, 0.76, 1.00, 0.00, "Beryllium"],
[  5, "B", 2.01, 0.84, 0.84, 1.92,  4,     10.811, 2.04,   8.298,   0.279723, 1.00, 0.71, 0.71, "Boron"],
[  6, "C", 2.50, 0.76, 0.76, 1.70,  4,    12.0107, 2.55, 11.2603,   1.262118, 0.40, 0.40, 0.40, "Carbon"],
[  7, "N", 3.07, 0.71, 0.71, 1.55,  4,    14.0067, 3.04, 14.5341,      -0.07, 0.05, 0.05, 1.00, "Nitrogen"],
[  8, "O", 3.50, 0.66, 0.66, 1.52,  2,    15.9994, 3.44, 13.6181,   1.461112, 1.00, 0.05, 0.05, "Oxygen"],
[  9, "F", 4.10, 0.57, 0.57, 1.47,  1, 18.9984032, 3.98, 17.4228,  3.4011887, 0.50, 0.70, 1.00, "Fluorine"],
[ 10, "Ne", 0.00, 0.58, 0.58, 1.54,  0,    20.1797, 0.00, 21.5645,          0, 0.70, 0.89, 0.96, "Neon"],
[ 11, "Na", 1.01, 1.66, 1.66, 2.27,  1,   22.98977, 0.93,  5.1391,   0.547926, 0.67, 0.36, 0.95, "Sodium"],
[ 12, "Mg", 1.23, 1.41, 1.41, 1.73,  2,    24.3050, 1.31,  7.6462,          0, 0.54, 1.00, 0.00, "Magnesium"],
[ 13, "Al", 1.47, 1.21, 1.21, 1.84,  6,  26.981538, 1.61,  5.9858,    0.43283, 0.75, 0.65, 0.65, "Aluminium"],
[ 14, "Si", 1.74, 1.11, 1.11, 2.10,  6,    28.0855, 1.90,  8.1517,   1.389521, 0.50, 0.60, 0.60, "Silicon"],
[ 15, "P", 2.06, 1.07, 1.07, 1.80,  6,  30.973761, 2.19, 10.4867,     0.7465, 1.00, 0.50, 0.00, "Phosphorus"],
[ 16, "S", 2.44, 1.05, 1.05, 1.80,  6,     32.065, 2.58,   10.36,  2.0771029, 0.70, 0.70, 0.00, "Sulfur"],
[ 17, "Cl", 2.83, 1.02, 1.02, 1.75,  1,     35.453, 3.16, 12.9676,   3.612724, 0.12, 0.94, 0.12, "Chlorine"],
[ 18, "Ar", 0.00, 1.06, 1.06, 1.88,  0,     39.948, 0.00, 15.7596,          0, 0.50, 0.82, 0.89, "Argon"],
[ 19, "K", 0.91, 2.03, 2.03, 2.75,  1,    39.0983, 0.82,  4.3407,   0.501459, 0.56, 0.25, 0.83, "Potassium"],
[ 20, "Ca", 1.04, 1.76, 1.76, 2.31,  2,     40.078, 1.00,  6.1132,    0.02455, 0.24, 1.00, 0.00, "Calcium"],
[ 21, "Sc", 1.20, 1.70, 1.70, 2.30,  6,   44.95591, 1.36,  6.5615,      0.188, 0.90, 0.90, 0.90, "Scandium"],
[ 22, "Ti", 1.32, 1.60, 1.60, 2.15,  6,     47.867, 1.54,  6.8281,      0.084, 0.75, 0.76, 0.78, "Titanium"],
[ 23, "V", 1.45, 1.53, 1.53, 2.05,  6,    50.9415, 1.63,  6.7462,      0.525, 0.65, 0.65, 0.67, "Vanadium"],
[ 24, "Cr", 1.56, 1.39, 1.39, 2.05,  6,    51.9961, 1.66,  6.7665,    0.67584, 0.54, 0.60, 0.78, "Chromium"],
[ 25, "Mn", 1.60, 1.39, 1.39, 2.05,  8,  54.938049, 1.55,   7.434,          0, 0.61, 0.48, 0.78, "Manganese"],
[ 26, "Fe", 1.64, 1.32, 1.32, 2.05,  6,     55.845, 1.83,  7.9024,      0.151, 0.88, 0.40, 0.20, "Iron"],
[ 27, "Co", 1.70, 1.26, 1.26, 2.00,  6,    58.9332, 1.88,   7.881,     0.6633, 0.94, 0.56, 0.63, "Cobalt"],
[ 28, "Ni", 1.75, 1.24, 1.24, 2.00,  6,    58.6934, 1.91,  7.6398,    1.15716, 0.31, 0.82, 0.31, "Nickel"],
[ 29, "Cu", 1.75, 1.32, 1.32, 2.00,  6,     63.546, 1.90,  7.7264,    1.23578, 0.78, 0.50, 0.20, "Copper"],
[ 30, "Zn", 1.66, 1.22, 1.22, 2.10,  6,      65.38, 1.65,  9.3942,          0, 0.49, 0.50, 0.69, "Zinc"],
[ 31, "Ga", 1.82, 1.22, 1.22, 1.87,  3,     69.723, 1.81,  5.9993,       0.41, 0.76, 0.56, 0.56, "Gallium"],
[ 32, "Ge", 2.02, 1.20, 1.20, 2.11,  4,      72.64, 2.01,  7.8994,   1.232712, 0.40, 0.56, 0.56, "Germanium"],
[ 33, "As", 2.20, 1.19, 1.19, 1.85,  3,   74.92160, 2.18,  9.7886,      0.814, 0.74, 0.50, 0.89, "Arsenic"],
[ 34, "Se", 2.48, 1.20, 1.20, 1.90,  2,      78.96, 2.55,  9.7524,    2.02067, 1.00, 0.63, 0.00, "Selenium"],
[ 35, "Br", 2.74, 1.20, 1.20, 1.83,  1,     79.904, 2.96, 11.8138,   3.363588, 0.65, 0.16, 0.16, "Bromine"],
[ 36, "Kr", 0.00, 1.16, 1.16, 2.02,  0,     83.798, 3.00, 13.9996,          0, 0.36, 0.72, 0.82, "Krypton"],
[ 37, "Rb", 0.89, 2.20, 2.20, 3.03,  1,    85.4678, 0.82,  4.1771,   0.485916, 0.44, 0.18, 0.69, "Rubidium"],
[ 38, "Sr", 0.99, 1.95, 1.95, 2.49,  2,      87.62, 0.95,  5.6949,    0.05206, 0.00, 1.00, 0.00, "Strontium"],
[ 39, "Y", 1.11, 1.90, 1.90, 2.40,  6,   88.90585, 1.22,  6.2173,      0.307, 0.58, 1.00, 1.00, "Yttrium"],
[ 40, "Zr", 1.22, 1.75, 1.75, 2.30,  6,     91.224, 1.33,  6.6339,      0.426, 0.58, 0.88, 0.88, "Zirconium"],
[ 41, "Nb", 1.23, 1.64, 1.64, 2.15,  6,   92.90638, 1.60,  6.7589,      0.893, 0.45, 0.76, 0.79, "Niobium"],
[ 42, "Mo", 1.30, 1.54, 1.54, 2.10,  6,      95.96, 2.16,  7.0924,     0.7472, 0.33, 0.71, 0.71, "Molybdenum"],
[ 43, "Tc", 1.36, 1.47, 1.47, 2.05,  6,         98, 1.90,    7.28,       0.55, 0.23, 0.62, 0.62, "Technetium"],
[ 44, "Ru", 1.42, 1.46, 1.46, 2.05,  6,     101.07, 2.20,  7.3605,    1.04638, 0.14, 0.56, 0.56, "Ruthenium"],
[ 45, "Rh", 1.45, 1.42, 1.42, 2.00,  6,  102.90550, 2.28,  7.4589,    1.14289, 0.04, 0.49, 0.55, "Rhodium"],
[ 46, "Pd", 1.35, 1.39, 1.39, 2.05,  6,     106.42, 2.20,  8.3369,    0.56214, 0.00, 0.41, 0.52, "Palladium"],
[ 47, "Ag", 1.42, 1.45, 1.45, 2.10,  6,   107.8682, 1.93,  7.5762,    1.30447, 0.88, 0.88, 1.00, "Silver"],
[ 48, "Cd", 1.46, 1.44, 1.44, 2.20,  6,    112.411, 1.69,  8.9938,          0, 1.00, 0.85, 0.56, "Cadmium"],
[ 49, "In", 1.49, 1.42, 1.42, 2.20,  3,    114.818, 1.78,  5.7864,      0.404, 0.65, 0.46, 0.45, "Indium"],
[ 50, "Sn", 1.72, 1.39, 1.39, 1.93,  4,    118.701, 1.96,  7.3439,   1.112066, 0.40, 0.50, 0.50, "Tin"],
[ 51, "Sb", 1.82, 1.39, 1.39, 2.17,  3,    121.760, 2.05,  8.6084,   1.047401, 0.62, 0.39, 0.71, "Antimony"],
[ 52, "Te", 2.01, 1.38, 1.38, 2.06,  2,     127.60, 2.10,  9.0096,   1.970875, 0.83, 0.48, 0.00, "Tellurium"],
[ 53, "I", 2.21, 1.39, 1.39, 1.98,  1,  126.90447, 2.66, 10.4513,   3.059038, 0.58, 0.00, 0.58, "Iodine"],
[ 54, "Xe", 0.00, 1.40, 1.40, 2.16,  0,    131.293, 2.60, 12.1298,          0, 0.26, 0.62, 0.69, "Xenon"],
[ 55, "Cs", 0.86, 2.44, 2.44, 3.43,  1,  132.90545, 0.79,  3.8939,   0.471626, 0.34, 0.09, 0.56, "Caesium"],
[ 56, "Ba", 0.97, 2.15, 2.15, 2.68,  2,    137.327, 0.89,  5.2117,    0.14462, 0.00, 0.79, 0.00, "Barium"],
[ 57, "La", 1.08, 2.07, 2.07, 2.50, 12,   138.9055, 1.10,  5.5769,       0.47, 0.44, 0.83, 1.00, "Lanthanum"],
[ 58, "Ce", 0.00, 2.04, 2.04, 2.48,  6,    140.116, 1.12,  5.5387,        0.5, 1.00, 1.00, 0.78, "Cerium"],
[ 59, "Pr", 0.00, 2.03, 2.03, 2.47,  6,  140.90765, 1.13,   5.473,        0.5, 0.85, 1.00, 0.78, "Praseodymium"],
[ 60, "Nd", 0.00, 2.01, 2.01, 2.45,  6,     144.24, 1.14,   5.525,        0.5, 0.78, 1.00, 0.78, "Neodymium"],
[ 61, "Pm", 0.00, 1.99, 1.99, 2.43,  6,        145, 0.00,   5.582,        0.5, 0.64, 1.00, 0.78, "Promethium"],
[ 62, "Sm", 0.00, 1.98, 1.98, 2.42,  6,     150.36, 1.17,  5.6437,        0.5, 0.56, 1.00, 0.78, "Samarium"],
[ 63, "Eu", 0.00, 1.98, 1.98, 2.40,  6,    151.964, 0.00,  5.6704,        0.5, 0.38, 1.00, 0.78, "Europium"],
[ 64, "Gd", 0.00, 1.96, 1.96, 2.38,  6,     157.25, 1.20,  6.1498,        0.5, 0.27, 1.00, 0.78, "Gadolinium"],
[ 65, "Tb", 0.00, 1.94, 1.94, 2.37,  6,  158.92534, 0.00,  5.8638,        0.5, 0.19, 1.00, 0.78, "Terbium"],
[ 66, "Dy", 0.00, 1.92, 1.92, 2.35,  6,    162.500, 1.22,  5.9389,        0.5, 0.12, 1.00, 0.78, "Dysprosium"],
[ 67, "Ho", 0.00, 1.92, 1.92, 2.33,  6,  164.93032, 1.23,  6.0215,        0.5, 0.00, 1.00, 0.61, "Holmium"],
[ 68, "Er", 0.00, 1.89, 1.89, 2.32,  6,    167.259, 1.24,  6.1077,        0.5, 0.00, 0.90, 0.46, "Erbium"],
[ 69, "Tm", 0.00, 1.90, 1.90, 2.30,  6,  168.93421, 1.25,  6.1843,        0.5, 0.00, 0.83, 0.32, "Thulium"],
[ 70, "Yb", 0.00, 1.87, 1.87, 2.28,  6,    173.054, 0.00,  6.2542,        0.5, 0.00, 0.75, 0.22, "Ytterbium"],
[ 71, "Lu", 0.00, 1.87, 1.87, 2.27,  6,   174.9668, 1.27,  5.4259,        0.5, 0.00, 0.67, 0.14, "Lutetium"],
[ 72, "Hf", 1.23, 1.75, 1.75, 2.25,  6,     178.49, 1.30,  6.8251,          0, 0.30, 0.76, 1.00, "Hafnium"],
[ 73, "Ta", 1.33, 1.70, 1.70, 2.20,  6,   180.9479, 1.50,  7.5496,      0.322, 0.30, 0.65, 1.00, "Tantalum"],
[ 74, "W", 1.40, 1.62, 1.62, 2.10,  6,     183.84, 2.36,   7.864,      0.815, 0.13, 0.58, 0.84, "Tungsten"],
[ 75, "Re", 1.46, 1.51, 1.51, 2.05,  6,    186.207, 1.90,  7.8335,       0.15, 0.15, 0.49, 0.67, "Rhenium"],
[ 76, "Os", 1.52, 1.44, 1.44, 2.00,  6,     190.23, 2.20,  8.4382,     1.0778, 0.15, 0.40, 0.59, "Osmium"],
[ 77, "Ir", 1.55, 1.41, 1.41, 2.00,  6,    192.217, 2.20,   8.967,    1.56436, 0.09, 0.33, 0.53, "Iridium"],
[ 78, "Pt", 1.44, 1.36, 1.36, 2.05,  6,    195.078, 2.28,  8.9588,     2.1251, 0.90, 0.85, 0.68, "Platinum"],
[ 79, "Au", 1.42, 1.36, 1.36, 2.10,  6,  196.96655, 2.54,  9.2255,    2.30861, 0.80, 0.82, 0.12, "Gold"],
[ 80, "Hg", 1.44, 1.32, 1.32, 2.05,  6,     200.59, 2.00, 10.4375,          0, 0.71, 0.71, 0.76, "Mercury"],
[ 81, "Tl", 1.44, 1.45, 1.45, 1.96,  3,   204.3833, 1.62,  6.1082,      0.377, 0.65, 0.33, 0.30, "Thallium"],
[ 82, "Pb", 1.55, 1.46, 1.46, 2.02,  4,      207.2, 2.33,  7.4167,      0.364, 0.34, 0.35, 0.38, "Lead"],
[ 83, "Bi", 1.67, 1.48, 1.48, 2.07,  3,  208.98040, 2.02,  7.2855,   0.942363, 0.62, 0.31, 0.71, "Bismuth"],
[ 84, "Po", 1.76, 1.40, 1.40, 1.97,  2,        209, 2.00,   8.414,        1.9, 0.67, 0.36, 0.00, "Polonium"],
[ 85, "At", 1.90, 1.50, 1.50, 2.02,  1,        210, 2.20,       0,        2.8, 0.46, 0.31, 0.27, "Astatine"],
[ 86, "Rn", 0.00, 1.50, 1.50, 2.20,  0,        222, 0.00, 10.7485,          0, 0.26, 0.51, 0.59, "Radon"],
[ 87, "Fr", 0.00, 2.60, 2.60, 3.48,  1,        223, 0.70,  4.0727,          0, 0.26, 0.00, 0.40, "Francium"],
[ 88, "Ra", 0.00, 2.21, 2.21, 2.83,  2,        226, 0.90,  5.2784,          0, 0.00, 0.49, 0.00, "Radium"],
[ 89, "Ac", 0.00, 2.15, 2.15, 2.00,  6,        227, 1.10,    5.17,          0, 0.44, 0.67, 0.98, "Actinium"],
[ 90, "Th", 0.00, 2.06, 2.06, 2.40,  6,   232.0381, 1.30,  6.3067,          0, 0.00, 0.73, 1.00, "Thorium"],
[ 91, "Pa", 0.00, 2.00, 2.00, 2.00,  6,  231.03588, 1.50,    5.89,          0, 0.00, 0.63, 1.00, "Protactinium"],
[ 92, "U", 0.00, 1.96, 1.96, 2.30,  6,  238.02891, 1.38,  6.1941,          0, 0.00, 0.56, 1.00, "Uranium"],
[ 93, "Np", 0.00, 1.90, 1.90, 2.00,  6,     237.05, 1.36,  6.2657,          0, 0.00, 0.50, 1.00, "Neptunium"],
[ 94, "Pu", 0.00, 1.87, 1.87, 2.00,  6,     244.06, 1.28,   6.026,          0, 0.00, 0.42, 1.00, "Plutonium"],
[ 95, "Am", 0.00, 1.80, 1.80, 2.00,  6,     243.06, 1.30,  5.9738,          0, 0.33, 0.36, 0.95, "Americium"],
[ 96, "Cm", 0.00, 1.69, 1.69, 2.00,  6,     247.07, 1.30,  5.9914,          0, 0.47, 0.36, 0.89, "Curium"],
[ 97, "Bk", 0.00, 1.60, 1.60, 2.00,  6,     247.07, 1.30,  6.1979,          0, 0.54, 0.31, 0.89, "Berkelium"],
[ 98, "Cf", 0.00, 1.60, 1.60, 2.00,  6,     251.08, 1.30,  6.2817,          0, 0.63, 0.21, 0.83, "Californium"],
[ 99, "Es", 0.00, 1.60, 1.60, 2.00,  6,     252.08, 1.30,    6.42,          0, 0.70, 0.12, 0.83, "Einsteinium"],
[100, "Fm", 0.00, 1.60, 1.60, 2.00,  6,     257.10, 1.30,     6.5,          0, 0.70, 0.12, 0.73, "Fermium"],
[101, "Md", 0.00, 1.60, 1.60, 2.00,  6,     258.10, 1.30,    6.58,          0, 0.70, 0.05, 0.65, "Mendelevium"],
[102, "No", 0.00, 1.60, 1.60, 2.00,  6,     259.10, 1.30,    6.65,          0, 0.74, 0.05, 0.53, "Nobelium"],
[103, "Lr", 0.00, 1.60, 1.60, 2.00,  6,     262.11, 0.00,     4.9,          0, 0.78, 0.00, 0.40, "Lawrencium"],
[104, "Rf", 0.00, 1.60, 1.60, 2.00,  6,     265.12, 0.00,       6,          0, 0.80, 0.00, 0.35, "Rutherfordium"],
[105, "Db", 0.00, 1.60, 1.60, 2.00,  6,     268.13, 0.00,       0,          0, 0.82, 0.00, 0.31, "Dubnium"],
[106, "Sg", 0.00, 1.60, 1.60, 2.00,  6,     271.13, 0.00,       0,          0, 0.85, 0.00, 0.27, "Seaborgium"],
[107, "Bh", 0.00, 1.60, 1.60, 2.00,  6,        270, 0.00,       0,          0, 0.88, 0.00, 0.22, "Bohrium"],
[108, "Hs", 0.00, 1.60, 1.60, 2.00,  6,     277.15, 0.00,       0,          0, 0.90, 0.00, 0.18, "Hassium"],
[109, "Mt", 0.00, 1.60, 1.60, 2.00,  6,     276.15, 0.00,       0,          0, 0.92, 0.00, 0.15, "Meitnerium"],
[110, "Ds", 0.00, 1.60, 1.60, 2.00,  6,     281.16, 0.00,       0,          0, 0.93, 0.00, 0.14, "Darmstadtium"],
[111, "Rg", 0.00, 1.60, 1.60, 2.00,  6,     280.16, 0.00,       0,          0, 0.94, 0.00, 0.13, "Roentgenium"],
[112, "Cn", 0.00, 1.60, 1.60, 2.00,  6,     285.17, 0.00,       0,          0, 0.95, 0.00, 0.12, "Copernicium"],
[113, "Nh", 0.00, 1.60, 1.60, 2.00,  6,     284.18, 0.00,       0,          0, 0.96, 0.00, 0.11, "Nihonium"],
[114, "Fl", 0.00, 1.60, 1.60, 2.00,  6,     289.19, 0.00,       0,          0, 0.97, 0.00, 0.10, "Flerovium"],
[115, "Mc", 0.00, 1.60, 1.60, 2.00,  6,     288.19, 0.00,       0,          0, 0.98, 0.00, 0.09, "Moscovium"],
[116, "Lv", 0.00, 1.60, 1.60, 2.00,  6,        293, 0.00,       0,          0, 0.99, 0.00, 0.08, "Livermorium"],
[117, "Ts", 0.00, 1.60, 1.60, 2.00,  6,        294, 0.00,       0,          0, 0.99, 0.00, 0.07, "Tennessine"],
[118, "Og", 0.00, 1.60, 1.60, 2.00,  6,        294, 0.00,       0,          0, 0.99, 0.00, 0.06, "Oganesson"],
]
element_list = []
'''Load the data from OpenBabel, and store it as both a
list of elements first, and then as an instance of Periodic Table.'''
for values in openbabel_element_data:
    number, symbol, AReneg, rcov, _, rvdw, maxbonds, MW, elneg, ionization, elaffinity, _, _, _, name = values
    AReneg = None if AReneg == 0 else AReneg
    rcov = None if rcov == 1.6 else rcov  # in Angstrom
    rvdw = None if rvdw == 2.0 else rvdw  # in Angstrom
    maxbonds = None if maxbonds == 6.0 else int(maxbonds)
    elneg = None if elneg == 0.0 else elneg
    ionization = None if ionization == 0.0 else ionization  # in eV
    elaffinity = None if elaffinity == 0.0 else elaffinity  # in eV
    index = number-1
    period = periods[index]
    group = groups[index]
    InChI_key = InChI_keys[index]
    cid = cids[index]
    phase = phases[index]
    Hf = Hfs[index]
    S0 = S0s[index]
    CAS = CAS_by_number[index]
    ele = Element(number, symbol, name, MW, CAS, AReneg, rcov, rvdw,
                  maxbonds, elneg, ionization, elaffinity, period, group,
                  cid, phase, Hf, S0, InChI_key)
    element_list.append(ele)

periodic_table = PeriodicTable(element_list)
'''Single instance of the PeriodicTable class. Use this, not the PeriodicTable
class directly.

A brief overview of using the periodic table and its elements:

>>> periodic_table.Na
<Element Sodium (Na) number 11, MW=22.98977>
>>> periodic_table.U.MW
238.02891
>>> periodic_table['Th'].CAS
'7440-29-1'
>>> periodic_table.lead.protons
82
>>> periodic_table['7440-57-5'].symbol
'Au'
>> len(periodic_table)
118
>>> 'gold' in periodic_table
True
>>> periodic_table.He.electrons, periodic_table.He.phase, periodic_table.He.Hf, periodic_table.He.S0
(2, 'g', 0.0, 126.2)
>>> periodic_table.Kr.block, periodic_table.Kr.period, periodic_table.Kr.group
('p', 4, 18)

'''
del openbabel_element_data

def molecular_weight(atoms):
    r'''Calculates molecular weight of a molecule given a dictionary of its
    atoms and their counts, in the format {symbol: count}.

    .. math::
        MW = \sum_i n_i MW_i

    Parameters
    ----------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]

    Returns
    -------
    MW : float
        Calculated molecular weight [g/mol]

    Notes
    -----
    Elemental data is from rdkit, with CAS numbers added. An exception is
    raised if an incorrect element symbol is given. Elements up to 118 are
    supported, as are deutreium and tritium.

    Examples
    --------
    >>> molecular_weight({'H': 12, 'C': 20, 'O': 5}) # DNA
    332.30628

    References
    ----------
    .. [1] RDKit: Open-source cheminformatics; http://www.rdkit.org
    
    '''
    MW = 0
    for i in atoms:
        if i in periodic_table:
            MW += periodic_table[i].MW*atoms[i]
        elif i == 'D':
            # Hardcoded MW until an actual isotope db is created
            MW += 2.014102*atoms[i]
        elif i == 'T':
            # Hardcoded MW until an actual isotope db is created
            MW += 3.0160492*atoms[i]
        else:
            raise ValueError('Molecule includes unknown atoms')
    return MW


def mass_fractions(atoms, MW=None):
    r'''Calculates the mass fractions of each element in a compound,
    given a dictionary of its atoms and their counts, in the format
    {symbol: count}.

    .. math::
        w_i =  \frac{n_i MW_i}{\sum_i n_i MW_i}

    Parameters
    ----------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]
    MW : float, optional
        Molecular weight, [g/mol]

    Returns
    -------
    mfracs : dict
        dictionary of mass fractions of individual atoms, indexed by symbol
        with proper capitalization, [-]

    Notes
    -----
    Molecular weight is optional, but speeds up the calculation slightly. It
    is calculated using the function `molecular_weight` if not specified.

    Elemental data is from rdkit, with CAS numbers added. An exception is
    raised if an incorrect element symbol is given. Elements up to 118 are
    supported.

    Examples
    --------
    >>> mass_fractions({'H': 12, 'C': 20, 'O': 5})
    {'H': 0.03639798802478244, 'C': 0.7228692758981262, 'O': 0.24073273607709128}

    References
    ----------
    .. [1] RDKit: Open-source cheminformatics; http://www.rdkit.org
    '''
    if not MW:
        MW = molecular_weight(atoms)
    mfracs = {}
    for i in atoms:
        if i in periodic_table:
            mfracs[i] = periodic_table[i].MW*atoms[i]/MW
        else:
            raise ValueError('Molecule includes unknown atoms')
    return mfracs


def atom_fractions(atoms):
    r'''Calculates the atomic fractions of each element in a compound,
    given a dictionary of its atoms and their counts, in the format
    {symbol: count}.

    .. math::
        a_i =  \frac{n_i}{\sum_i n_i}

    Parameters
    ----------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]

    Returns
    -------
    afracs : dict
        dictionary of atomic fractions of individual atoms, indexed by symbol
        with proper capitalization, [-]

    Notes
    -----
    No actual data on the elements is used, so incorrect or custom compounds
    would not raise an error.

    Examples
    --------
    >>> atom_fractions({'H': 12, 'C': 20, 'O': 5})
    {'H': 0.32432432432432434, 'C': 0.5405405405405406, 'O': 0.13513513513513514}

    References
    ----------
    .. [1] RDKit: Open-source cheminformatics; http://www.rdkit.org
    '''
    count = sum(atoms.values())
    afracs = {}
    for i in atoms:
        afracs[i] = atoms[i]/count
    return afracs


def mixture_atomic_composition(atomss, zs):
    r'''Simple function to calculate the atomic average composition of a
    mixture, using the mole fractions of each species and their own atomic 
    compositions.
    
    Parameters
    ----------
    atomss : list[dict[(str, int)]]
        List of dictionaries of atomic compositions, [-]
    zs : list[float]
        Mole fractions of each component, [-]

    Returns
    -------
    atoms : dict[(str, int)]
        Atomic composition

    Notes
    -----

    Examples
    --------
    >>> mixture_atomic_composition([{'O': 2}, {'N': 1, 'O': 2}, {'C': 1, 'H': 4}], [0.95, 0.025, .025])
    {'O': 1.95, 'N': 0.025, 'C': 0.025, 'H': 0.1}
    '''
    ans = {}
    for atoms, zs_i in zip(atomss, zs):
        for key, val in atoms.items():
            if key in ans:
                ans[key] += val*zs_i
            else:
                ans[key] = val*zs_i
    return ans


def mixture_atomic_composition_ordered(atomss, zs):
    r'''Simple function to calculate the atomic average composition of a
    mixture, using the mole fractions of each species and their own atomic 
    compositions. Returns the result as a sorted list with atomic numbers from
    low to high.
    
    Parameters
    ----------
    atomss : list[dict[(str, int)]]
        List of dictionaries of atomic compositions, [-]
    zs : list[float]
        Mole fractions of each component; this can also be a molar flow rate
        and then the `abundances` will be flows, [-]

    Returns
    -------
    abundances : list[float]
        Number of atoms of each element per mole of the feed, [-]
    atom_keys : list[str]
        Atomic elements, sorted from lowest atomic number to highest

    Notes
    -----
    Useful to ensure a matrix order is consistent in multiple steps.

    Examples
    --------
    >>> mixture_atomic_composition_ordered([{'O': 2}, {'N': 1, 'O': 2}, {'C': 1, 'H': 4}], [0.95, 0.025, .025])
    ([0.1, 0.025, 0.025, 1.95], ['H', 'C', 'N', 'O'])
    '''
    ans = mixture_atomic_composition(atomss, zs)
    nums = []
    eles = []
    for k, n in sorted(ans.items(), key=lambda x: periodic_table[x[0]].number):
        nums.append(n)
        eles.append(k)
    return nums, eles


def atom_matrix(atomss, atom_IDs=None):
    r'''Simple function to create a matrix of elements in each compound, where
    each row has the same elements.
    
    Parameters
    ----------
    atomss : list[dict[(str, int)]]
        List of dictionaries of atomic compositions, [-]
    atom_IDs : list[str], optional
        Optionally, a subset (or simply ordered differently) of elements to 
        consider, [-]

    Returns
    -------
    matrix : list[list[float]]
        The number of each element in each compound as a matrix, indexed as
        [compound][element], [-]

    Notes
    -----

    Examples
    --------
    >>> atom_matrix([{'C': 1, 'H': 4}, {'C': 2, 'H': 6}, {'N': 2}, {'O': 2}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}])
    [[4, 1, 0.0, 0.0], [6, 2, 0.0, 0.0], [0.0, 0.0, 2, 0.0], [0.0, 0.0, 0.0, 2], [2, 0.0, 0.0, 1], [0.0, 1, 0.0, 2]]
    '''
    if atom_IDs is None:
        ans = set([])
        for atoms in atomss:
            for i in atoms.keys():
                ans.add(i)
        atom_IDs = sorted(list(ans), key=lambda x: periodic_table[x].number)
    
    atom_idx = {k: i for i, k in enumerate(atom_IDs)}
    n_atoms = len(atom_IDs)
    element_matrix = []
    for atoms in atomss:
        l = [0.0]*n_atoms
        for k, v in atoms.items():
            try:
                l[atom_idx[k]] = v
            except KeyError:
                pass
        element_matrix.append(l)
        
    return element_matrix

def similarity_variable(atoms, MW=None):
    r'''Calculates the similarity variable of an compound, as defined in [1]_.
    Currently only applied for certain heat capacity estimation routines.

    .. math::
        \alpha = \frac{N}{MW} = \frac{\sum_i n_i}{\sum_i n_i MW_i}

    Parameters
    ----------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]
    MW : float, optional
        Molecular weight, [g/mol]

    Returns
    -------
    similarity_variable : float
        Similarity variable as defined in [1]_, [mol/g]

    Notes
    -----
    Molecular weight is optional, but speeds up the calculation slightly. It
    is calculated using the function `molecular_weight` if not specified.

    Examples
    --------
    >>> similarity_variable({'H': 32, 'C': 15})
    0.2212654140784498

    References
    ----------
    .. [1] Laštovka, Václav, Nasser Sallamie, and John M. Shaw. "A Similarity
       Variable for Estimating the Heat Capacity of Solid Organic Compounds:
       Part I. Fundamentals." Fluid Phase Equilibria 268, no. 1-2
       (June 25, 2008): 51-60. doi:10.1016/j.fluid.2008.03.019.
    '''
    if not MW:
        MW = molecular_weight(atoms)
    return sum(atoms.values())/MW

def atoms_to_Hill(atoms):
    r'''Determine the Hill formula of a compound, given a dictionary of its
    atoms and their counts, in the format {symbol: count}.

    Parameters
    ----------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]

    Returns
    -------
    Hill_formula : str
        Hill formula, [-]

    Notes
    -----
    The Hill system is as follows:

    If the chemical has 'C' in it, this is listed first, and then if it has
    'H' in it as well as 'C', then that goes next. All elements are sorted
    alphabetically afterwards, including 'H' if 'C' is not present.
    All elements are followed by their count, unless it is 1.

    Examples
    --------
    >>> atoms_to_Hill({'H': 5, 'C': 2, 'Br': 1})
    'C2H5Br'

    References
    ----------
    .. [1] Hill, Edwin A."“ON A SYSTEM OF INDEXING CHEMICAL LITERATURE;
       ADOPTED BY THE CLASSIFICATION DIVISION OF THE U. S. PATENT OFFICE.1."
       Journal of the American Chemical Society 22, no. 8 (August 1, 1900):
       478-94. doi:10.1021/ja02046a005.
    '''
    def str_ele_count(ele):
        if atoms[ele] == 1:
            count = ''
        else:
            count = str(atoms[ele])
        return count
    atoms = atoms.copy()
    s = ''
    if 'C' in atoms.keys():
        s += 'C' + str_ele_count('C')
        del atoms['C']
        if 'H' in atoms.keys():
            s += 'H' + str_ele_count('H')
            del atoms['H']
        for ele in sorted(atoms.keys()):
            s += ele + str_ele_count(ele)
    else:
        for ele in sorted(atoms.keys()):
            s += ele + str_ele_count(ele)
    return s

_simple_formula_parser_re_str = r'([A-Z][a-z]{0,2})([\d\.\d]+)?'
_simple_formula_parser_re = None # Delay creation to simple_formula_parser to speedup start

def simple_formula_parser(formula):
    r'''Basic formula parser, primarily for obtaining element counts from 
    formulas as formated in PubChem. Handles formulas with integer or decimal
    counts (with period separator), but no brackets, no hydrates, no charges,
    no isotopes, and no group multipliers.
    
    Strips charges from the end of a formula first. Accepts repeated chemical
    units. Performs no sanity checking that elements are actually elements.
    As it uses regular expressions for matching, errors are mostly just ignored.
    
    Parameters
    ----------
    formula : str
        Formula string, very simply formats only.

    Returns
    -------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]

    Notes
    -----
    Inspiration taken from the thermopyl project, at
    https://github.com/choderalab/thermopyl.

    Examples
    --------
    >>> simple_formula_parser('CO2')
    {'C': 1, 'O': 2}
    '''
    global _simple_formula_parser_re
    if not _simple_formula_parser_re:
        _simple_formula_parser_re = re.compile(_simple_formula_parser_re_str)
    formula = formula.split('+')[0].split('-')[0]
    counts = {}
    for element, count in _simple_formula_parser_re.findall(formula):
        if count.isdigit():
            count = int(count)
        elif count:
            count = float(count)
        else:
            count = 1
        if element in counts:
            counts[element] += count
        else:
            counts[element] = count
    return counts

#  Delay creation to simple_formula_parser to speedup start
formula_token_matcher_rational_re_str = r'[A-Z][a-z]?|(?:\d*[.])?\d+|\d+|[()]'
bracketed_charge_re_str = r'\([+-]?\d+\)$|\(\d+[+-]?\)$|\([+-]+\)$'
formula_token_matcher_rational = bracketed_charge_re = None
letter_set = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')

def nested_formula_parser(formula, check=True):
    r'''Improved formula parser which handles braces and their multipliers, 
    as well as rational element counts.

    Strips charges from the end of a formula first. Accepts repeated chemical
    units. Performs no sanity checking that elements are actually elements.
    As it uses regular expressions for matching, errors are mostly just ignored.
    
    Parameters
    ----------
    formula : str
        Formula string, very simply formats only.
    check : bool
        If `check` is True, a simple check will be performed to determine if
        a formula is not a formula and an exception will be raised if it is
        not, [-]

    Returns
    -------
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]

    Notes
    -----
    Inspired by the approach taken by CrazyMerlyn on a reddit DailyProgrammer
    challenge, at https://www.reddit.com/r/dailyprogrammer/comments/6eerfk/20170531_challenge_317_intermediate_counting/

    Examples
    --------
    >>> nested_formula_parser('Pd(NH3)4.0001+2')
    {'Pd': 1, 'N': 4.0001, 'H': 12.0003}
    '''
    global formula_token_matcher_rational, bracketed_charge_re
    if formula_token_matcher_rational is None:
        formula_token_matcher_rational = re.compile(formula_token_matcher_rational_re_str)
        bracketed_charge_re = re.compile(bracketed_charge_re_str)
    
    formula = formula.replace('[', '').replace(']', '')
    charge_splits = bracketed_charge_re.split(formula)
    if len(charge_splits) > 1:
        formula = charge_splits[0]
    else:
        formula = formula.split('+')[0].split('-')[0]
    
    stack = [[]]
    last = stack[0]
    tokens = formula_token_matcher_rational.findall(formula)
    # The set of letters in the tokens should match the set of letters
    if check:
        token_letters = set([j for i in tokens for j in i if j in letter_set])
        formula_letters = set(i for i in formula if i in letter_set)
        if formula_letters != token_letters:
            raise ValueError('Input may not be a formula; extra letters were detected')
    
    for token in tokens:
        if token == "(":
            stack.append([])
            last = stack[-1]
        elif token == ")":
            temp_dict = {}
            for d in last:
                for ele, count in d.items():
                    if ele in temp_dict:
                        temp_dict[ele] = temp_dict[ele] + count
                    else:
                        temp_dict[ele] = count
            stack.pop()
            last = stack[-1]
            last.append(temp_dict)
        elif token.isalpha():
            last.append({token: 1})
        else:
            v = float(token)
            v_int = int(v)
            if v_int == v:
                v = v_int
            last[-1] = {ele: count*v for ele, count in last[-1].items()}
    ans = {}
    for d in last:
        for ele, count in d.items():
            if ele in ans:
                ans[ele] = ans[ele] + count
            else:
                ans[ele] = count
    return ans

def charge_from_formula(formula):
    r'''Basic formula parser to determine the charge from a formula - given
    that the charge is already specified as one element of the formula.

    Performs no sanity checking that elements are actually elements.
    
    Parameters
    ----------
    formula : str
        Formula string, very simply formats only, ending in one of '+x',
        '-x', n*'+', or n*'-' or any of them surrounded by brackets but always
        at the end of a formula.

    Returns
    -------
    charge : int
        Charge of the molecule, [faraday]

    Notes
    -----

    Examples
    --------
    >>> charge_from_formula('Br3-')
    -1
    >>> charge_from_formula('Br3(-)')
    -1
    '''
    global bracketed_charge_re
    negative = '-' in formula
    positive = '+' in formula
    if positive and negative:
        raise ValueError('Both negative and positive signs were found in the formula; only one sign is allowed')
    elif not (positive or negative):
        return 0
    multiplier, sign = (-1, '-') if negative else (1, '+')
    
    hit = False
    if '(' in formula:
        if bracketed_charge_re is None: # pragma: no cover
            bracketed_charge_re = re.compile(bracketed_charge_re_str)
        hit = bracketed_charge_re.findall(formula)
        if hit:
            formula = hit[-1].replace('(', '').replace(')', '')

    count = formula.count(sign)
    if count == 1:
        splits = formula.split(sign)
        if splits[1] == '' or splits[1] == ')':
            return multiplier
        else:
            return multiplier*int(splits[1])
    else:
        return multiplier*count

def serialize_formula(formula):
    r'''Basic formula serializer to construct a consistently-formatted formula.
    This is necessary for handling user-supplied formulas, which are not always
    well formatted.

    Performs no sanity checking that elements are actually elements.
    
    Parameters
    ----------
    formula : str
        Formula string as parseable by the method nested_formula_parser, [-]

    Returns
    -------
    formula : str
        A consistently formatted formula to describe a molecular formula, [-]

    Notes
    -----

    Examples
    --------
    >>> serialize_formula('Pd(NH3)4+3')
    'H12N4Pd+3'
    '''
    charge = charge_from_formula(formula)
    element_dict = nested_formula_parser(formula)
    base = atoms_to_Hill(element_dict)
    if charge == 0:
        pass
    elif charge > 0:
        if charge == 1:
            base += '+'
        else:
            base += '+' + str(charge)
    elif charge < 0:
        if charge == -1:
            base += '-'
        else:
            base +=  str(charge)
    return base
    
