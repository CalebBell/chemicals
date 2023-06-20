import json
import os

import pandas as pd

from chemicals.identifiers import check_CAS

dirname = os.path.dirname(os.path.abspath(__file__))

metadata_table = pd.read_csv(os.path.join(dirname, 'janaf_tables_index.csv'), index_col=0)
# print(metadata_table)
# Identify the CAS numbers, that, well, aren't
bad_CASs = set()
for CAS in metadata_table.index.values.tolist():
    try:
        assert check_CAS(CAS)
    except:
        bad_CASs.add(CAS)
        # print(CAS)

"""This section is from the thermochem project:
With changes to export the data.


Copyright (c) 2007-2008 by the respective authors (see AUTHORS file).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * The names of the contributors may not be used to endorse or
      promote products derived from this software without specific
      prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

"""
This module gets thermodynamic data from the JANAF database.
Files are downloaded from the NIST servers as needed and then cached locally.

Zack Gainsforth

Funding by NASA
"""

import os
import sys
from textwrap import dedent

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

try:
    # Python 3
    import urllib.request as urllib2
except ImportError:
    # Python 2
    import urllib2

try:
    # Python 2
    from StringIO import StringIO
except ImportError:
    # Python 3
    from io import StringIO

# Universal gas constant R
R = 8.314472
# Reference temp
Tr = 298.15 # K


class JanafPhase:
    """
    Class which is created by Janafdb for a specific phase.

    It reads in the JANAF data file and produces functions which interpolate
    the thermodynamic constants.

    Tr stands for reference temperature and is 298.15 K

    >>> db = Janafdb()
    >>> p = db.getphasedata(formula='O2Ti', name='Rutile', phase='cr')
    >>> p.cp([500, 550, 1800]).astype(int).tolist()
    [67, 68, 78]
    >>> print(p.S([500, 550, 1800]))        # Entropy in J/mol/K
    [  82.201    88.4565  176.876 ]
    >>> print(p.gef([500, 550, 1800]))      # [G-H(Tr)]/T in J/mol/K
    [  57.077   59.704  115.753]
    >>> print(p.hef([500, 550, 1800]))      # H-H(Tr) in kJ/mol
    [  12.562    15.9955  110.022 ]
    >>> print(p.DeltaH([500, 550, 1800]))   # Standard enthalpy of formation in kJ/mol
    [-943.670  -943.2295 -936.679 ]
    >>> print(p.DeltaG([500, 550, 1800]))   # Gibbs free enegy in kJ/mol
    [-852.157  -843.0465 -621.013 ]
    >>> p.logKf([500, 550, 1800]).astype(int).tolist() # Equilibrium constant of formation.
    [89, 80, 18]
    >>> print(p.cp(1000))                   # Heat capacity in J/mol/K
    74.852
    >>> print(p.cp(50000))                  # Example of erroneous extrapolation.
    Traceback (most recent call last):
        ...
    ValueError: A value in x_new is above the interpolation range.
    """

    def __init__(self, rawdata_text):
        # Store the raw data text file from NIST.
        self.rawdata_text = rawdata_text

        self.description = self.rawdata_text.splitlines()[0]

        # Read the text file into a DataFrame.
        # TODO: adjust usecols length to be within bounds, Pandas deprecation
        data = pd.read_csv(
            StringIO(self.rawdata_text),
            skiprows=2,
            header=None,
            delimiter=r'[\t\s]+',
            engine='python',
            names=['T', 'Cp', 'S', '[G-H(Tr)]/T', 'H-H(Tr)', 'Delta_fH', 'Delta_fG', 'log(Kf)'],
            usecols=['T', 'Cp', 'S', '[G-H(Tr)]/T', 'H-H(Tr)', 'Delta_fH', 'Delta_fG', 'log(Kf)'] # Ignore extra columns -- those are caused by comments in the text file
        )
        self.rawdata = data

        # Sometimes the JANAF files have funky stuff written in them.
        # (Old school text format...)
        # Clean it up.
        for c in data.columns:
            # We only need to polish up columns that aren't floating point
            # numbers.
            if np.issubdtype(data.dtypes[c], np.floating):
                continue
            # Change INFINITE to inf
            data.loc[data[c] == 'INFINITE', c] = np.inf
            # Anything else becomes a nan.
            # Convert to floats.
            data[c] = pd.to_numeric(data[c], errors='coerce')

        # Handle NaNs for the phase transition points. This only affects
        # Delta_fG, Delta_fH, and log(Kf)
        good_indices = np.where(np.isfinite(data['Delta_fH']))

        # Now make interpolatable functions for each of these.
        self.cp = interp1d(self.rawdata['T'], self.rawdata['Cp'])
        self.S = interp1d(self.rawdata['T'], self.rawdata['S'])
        self.gef = interp1d(self.rawdata['T'], self.rawdata['[G-H(Tr)]/T'])
        self.hef = interp1d(self.rawdata['T'], self.rawdata['H-H(Tr)'])
        self.DeltaH = interp1d(self.rawdata['T'].iloc[good_indices],
                               self.rawdata['Delta_fH'].iloc[good_indices])
        self.DeltaG = interp1d(self.rawdata['T'].iloc[good_indices],
                               self.rawdata['Delta_fG'].iloc[good_indices])
        self.logKf = interp1d(self.rawdata['T'].iloc[good_indices],
                              self.rawdata['log(Kf)'].iloc[good_indices])

    def __str__(self):
        rep = super().__str__()
        rep += "\n  "
        rep += self.description
        rep += f"\n    Cp({Tr:0.2f}) = {self.cp(Tr):0.3f} J/mol/K"
        rep += f"\n    S({Tr:0.2f}) = {self.S(Tr):0.3f} J/mol/K"
        rep += f"\n    [G-H({Tr:0.2f})]/{Tr:0.2f} = {self.gef(Tr):0.3f} J/mol/K"
        rep += f"\n    H-H({Tr:0.2f}) = {self.hef(Tr):0.3f} J/mol/K"
        rep += f"\n    Delta_fH({Tr:0.2f}) = {self.DeltaH(Tr):0.3f} kJ/mol"
        rep += f"\n    Delta_fG({Tr:0.2f}) = {self.DeltaG(Tr):0.3f} kJ/mol"
        rep += f"\n    log(Kf(({Tr:0.2f})) = {self.logKf(Tr):0.3f}"
        return rep


class Janafdb:
    """
    Class that reads the NIST JANAF tables for thermodynamic data.

    Data is initially read from the web servers, and then cached.

    Examples
    --------
    >>> rutile = Janafdb().getphasedata(name='Rutile')

    To load thermodynamic constants for TiO2, rutile.
    """

    VALIDPHASETYPES = ['cr', 'l', 'cr,l', 'g', 'ref', 'cd', 'fl', 'am', 'vit',
                       'mon', 'pol', 'sln', 'aq', 'sat']
    JANAF_URL = "https://janaf.nist.gov/tables/%s.txt"

    def __init__(self):
        """
        We have an index file which can be used to build the url for all phases
        on the NIST site.
        """
        # Read the index file which tells us the filenames for all the phases
        # in the JANAF database.
        __file__
        janaf_index =os.path.join(dirname, 'JANAF_index.txt')
        self.db = pd.read_csv(janaf_index, delimiter='|', header=None)
        # Name the columns and trim whitespace off the text fields.
        self.db.columns = ['formula', 'name', 'phase', 'filename']
        self.db["formula"] = self.db["formula"].map(str.strip)
        self.db["name"] = self.db["name"].map(str.strip)
        self.db["phase"] = self.db["phase"].map(str.strip)
        self.db["filename"] = self.db["filename"].map(str.strip)

        # Make sure that the directory for cached JANAF files exists.
        self.JANAF_cachedir = os.path.join(dirname, 'JANAF_Cache')
        if not os.path.exists(self.JANAF_cachedir):
            os.mkdir(self.JANAF_cachedir)

    def __str__(self):
        rep = super().__str__()
        # rep = "\tFormula = %s"%self.db["formula"]
        rep += "\n  Try:\n"
        rep += "    Janafdb().search('Ti')\n"
        rep += "  or:\n"
        rep += "    Janafdb().getphasedata(name='Magnesium Oxide', phase='l')\n"
        return rep

    def search(self, searchstr):
        """
        List all the species containing a string. Helpful for
        interactive use of the database.

        Parameters
        ----------
        searchstr : str
            The search string to look for

        Returns
        -------
        pandas.DataFrame
            Dataframe containing valid phases

        Examples
        --------
        >>> db = Janafdb()
        >>> s = db.search('Rb-')
        >>> print(s)
             formula           name phase filename
        1710     Rb-  Rubidium, Ion     g   Rb-007
        >>> s = db.search('Ti')
        >>> print(len(s))
        88
        """
        formulasearch = self.db['formula'].str.contains(searchstr)
        namesearch = self.db['name'].str.contains(searchstr)

        return self.db[formulasearch | namesearch]

    def getphasedata(self,
                     formula=None,
                     name=None,
                     phase=None,
                     filename=None,
                     cache=True):
        """
        Returns an element instance given the name of the element.
        formula, name and phase match the respective fields in the JANAF index.

        Parameters
        ----------
        formula : str
            Select records that match the chemical formula
        name : str
            Select records that match the chemical/mineral name
        phase : str
            Select records that match the chemical phase.
            Must be one of the following valid phases:
            cr, l, cr,l, g, ref, cd, fl, am, vit, mon, pol, sln, aq, sat
        filename : str
            Select only records that match the filename on the website, which
            is very unique.
        cache : bool, default True
            Whether to cache the Janaf database. Setting this to false will
            download the Janaf database every time it is used.

        Examples
        --------
        >>> db = Janafdb()
        >>> db.getphasedata(formula='O2Ti', phase='cr')
        Traceback (most recent call last):
            ...
        ValueError: There are 2 records matching this pattern:
            ...
        Please select a unique record.
        >>> db.getphasedata(formula='Oxyz')
        Traceback (most recent call last):
            ...
        ValueError: Did not find a phase with formula = Oxyz
                    Please provide enough information to select a unique record.
                    Also check that you didn't eliminate the record you want by choosing too many constraints where one or more constraint is incorrect.
        >>> db.getphasedata(formula='Oxyz', phase='l')
        Traceback (most recent call last):
            ...
        ValueError: Did not find a phase with formula = Oxyz, phase = l
                    Please provide enough information to select a unique record.
                    Also check that you didn't eliminate the record you want by choosing too many constraints where one or more constraint is incorrect.
        >>> FeO = db.getphasedata(formula='FeO', phase='cr,l')
        >>> print(FeO)
        <thermochem.janaf.JanafPhase object at 0x...>
          Iron Oxide (FeO)  Fe1O1(cr,l)
            Cp(298.15) = 49.915 J/mol/K
            S(298.15) = 60.752 J/mol/K
            [G-H(298.15)]/298.15 = 60.752 J/mol/K
            H-H(298.15) = 0.000 J/mol/K
            Delta_fH(298.15) = -272.044 kJ/mol
            Delta_fG(298.15) = -251.429 kJ/mol
            log(Kf((298.15)) = 44.049
        """
        # Check that the phase type requested is valid.
        if phase is not None:
            phase = phase.lower()
        if phase is not None and phase not in self.VALIDPHASETYPES:
            raise ValueError("Valid phase types are %s." % self.VALIDPHASETYPES)

        # We can search on either an exact formula, partial text match in the
        # name, and exact phase type.
        # Start with all records selected in the search, and then we trim.
        formulasearch = pd.Series(np.ones(len(self.db)), dtype=bool)
        namesearch = formulasearch.copy()
        phasesearch = formulasearch.copy()
        filenamesearch = formulasearch.copy()
        if formula is not None:
            formulasearch = self.db['formula'] == formula
        if name is not None:
            namesearch = self.db['name'].str.lower() == name.lower()
            # namesearch = self.db['name'].str.lower().str.equals(name.lower(), regex=False)
        if phase is not None:
            phasesearch = self.db['phase'] == phase
        if filename is not None:
            filenamesearch = self.db['filename'].str.lower() == filename.lower()
        # Combine.
        searchmatch = formulasearch & namesearch & phasesearch & filenamesearch

        # Get the record (should be one record) which specifies this phase.
        phase_record = self.db[searchmatch]
        if phase_record.empty:
            searched = []
            if formula is not None:
                searched.append("formula = %s" % formula)
            if phase is not None:
                searched.append("phase = %s" % phase)
            if filename is not None:
                searched.append("filename = %s" % filename)
            search_string = ", ".join(searched)
            raise ValueError("""
            Did not find a phase with %s
            Please provide enough information to select a unique record.
            Also check that you didn't eliminate the record you want by choosing too many constraints where one or more constraint is incorrect.""" % search_string)
        if len(phase_record) > 1:
            # The user has entered in data that does not uniquely select one
            # record. Let's help him out by listing his options unless it is
            # too many.
            raise ValueError(dedent("""
            There are %d records matching this pattern:
            %s

            Please select a unique record.""") % (len(phase_record), phase_record))

        # At this point we have one record.  Check if we have that file cached.
        cachedfilename = os.path.join(
            self.JANAF_cachedir,
            "%s.txt" % phase_record['filename'].values[0]
        )
        if cache and os.path.exists(cachedfilename):
            # Yes it was cached, so let's read it into memory.
            with open(cachedfilename) as f:
                textdata = f.read()
        else:
            # No it was not cached so let's get it from the web.
            response = urllib2.urlopen(Janafdb.JANAF_URL %
                                       phase_record['filename'].values[0])
            textdata = response.read()
            if sys.version_info[0] > 2:
                textdata = textdata.decode()

            # And cache the data so we aren't making unnecessary trips to the
            # web.
            if cache:
                with open(cachedfilename, 'w') as f:
                    f.write(textdata)

        # Create a phase class and return it.
        return JanafPhase(textdata)


db = Janafdb()

# Get as much data our for gases and liquids as possible (easy to understand)
# Not copyright as above
from math import isnan

Hfgs_dict = {}
S0gs_dict = {}
Gfgs_dict = {}
Cpgs_dict = {}
Cpgs_values_dict = {}
Tg_values_dict = {}
names_dict = {}
formulas_dict = {}

Hfls_dict = {}
S0ls_dict = {}
Gfls_dict = {}
Cpls_dict = {}
Cpls_values_dict = {}
Tl_values_dict = {}

Hfss_dict = {}
S0ss_dict = {}
Gfss_dict = {}
Cpss_dict = {}
Cpss_values_dict = {}
Ts_values_dict = {}

bad_states = {'l,g', 'cr,l'}

solid_CAS_trans = {
# (CAS, name, state) -> Fake CAS
('7440-67-7', 'Zirconium, Alpha', 'cr'): '2099592000-00-0',
('7440-67-7', 'Zirconium, Beta', 'cr'): '2099576000-00-0',
# ('7782-42-5', 'Carbon', 'ref'): '7782-42-5',

}

from collections import Counter

cr_CAS_counter = Counter()
for CAS, row in metadata_table.iterrows():
    formula = row['Formula']
    name = row['Name']
    state = row['State']
    if state == 'cr':
        cr_CAS_counter[CAS] += 1


for CAS, row in metadata_table.iterrows():
    if CAS in bad_CASs:
        continue
    # print('Processing CAS', CAS)
    formula = row['Formula']
    name = row['Name']
    state = row['State']
    is_cr_unique = cr_CAS_counter[CAS] == 1
    is_graphite = CAS == '7440-44-0' and state == 'ref'

    solid_key = (CAS, name, state)

    if solid_key in solid_CAS_trans:
        CAS = solid_CAS_trans[solid_key]
        # print('changed!')

    if state in bad_states and not is_graphite:
        continue

    if name == 'Nitrous Acid, Cis':
        # Same CAS as Nitrous Acid, Trans but that one is more stable
        continue
    if CAS == '7664-39-3' and formula != 'HF':
        # H2F2 to H7F7 as well not sure why they are there
        continue
    if name == 'Iron Carbonyl':
        # Incorrectly has liquid and solid data, cannot use it
        continue
    if ' bar' in name.lower():
        # bad waters
        continue


    names_dict[CAS] = name
    formulas_dict[CAS] = formula

    if state not in ('l', 'g'):
        pass
        # print(CAS, formula, name, state)
    if state == 'cr' and not is_cr_unique :
        print([CAS, formula, name, state])

    try:
        p = db.getphasedata(formula=formula, name=name, phase=state)
    except Exception as e:
        print(['Failed', CAS, formula, name, state], e)
        p = db.getphasedata(formula=formula, name=name, phase=state)
        continue
    Hf_298 = float(p.DeltaH([298.15])[0])*1000.0
    Gf_298 = float(p.DeltaG([298.15])[0])*1000.0
    S0_298 =  float(p.S([298.15])[0])
    Cp_298 =  float(p.cp([298.15])[0])
    Ts = p.rawdata['T'].tolist()
    Cps = p.rawdata['Cp'].tolist()
    Ts_Cp = []
    Cps_Cps = []
    for T, Cpi in zip(Ts, Cps):
        if not isnan(T) and not isnan(Cpi):
            Ts_Cp.append(T)
            Cps_Cps.append(Cpi)
    Ts = Ts_Cp
    Cps = Cps_Cps

    assert not isnan(Hf_298)
    assert not isnan(Gf_298)
    assert not isnan(S0_298)
    assert not isnan(Cp_298)
    if state == 'g':
        Hfgs_dict[CAS] = Hf_298
        Gfgs_dict[CAS] = Gf_298
        S0gs_dict[CAS] = S0_298
        Cpgs_dict[CAS] = Cp_298
        assert CAS not in Cpgs_values_dict, CAS
        if len(Ts_Cp) == len(set(Ts_Cp)) or True:
            # Skip duplicate T values as the table isn't really single phase
            Cpgs_values_dict[CAS] = Cps
            Tg_values_dict[CAS] = Ts_Cp
    elif state == 'l':
        Hfls_dict[CAS] = Hf_298
        Gfls_dict[CAS] = Gf_298
        S0ls_dict[CAS] = S0_298
        Cpls_dict[CAS] = Cp_298
        assert CAS not in Cpls_values_dict, CAS
        if len(Ts_Cp) == len(set(Ts_Cp)) or True:
            # Skip duplicate T values as the table isn't really single phase
            Cpls_values_dict[CAS] = Cps
            Tl_values_dict[CAS] = Ts_Cp
    elif state == 'cr' and (solid_key in solid_CAS_trans or is_cr_unique) or is_graphite:
        Hfss_dict[CAS] = Hf_298
        Gfss_dict[CAS] = Gf_298
        S0ss_dict[CAS] = S0_298
        Cpss_dict[CAS] = Cp_298
        assert CAS not in Cpss_values_dict, CAS
        if len(Ts_Cp) == len(set(Ts_Cp)) or True:
            # Skip duplicate T values as the table isn't really single phase
            Cpss_values_dict[CAS] = Cps
            Ts_values_dict[CAS] = Ts_Cp

        if is_graphite:
            graphite_CAS = '7782-42-5'
            Hfss_dict[graphite_CAS] = Hf_298
            Gfss_dict[graphite_CAS] = Gf_298
            S0ss_dict[graphite_CAS] = S0_298
            Cpss_dict[graphite_CAS] = Cp_298
            Cpss_values_dict[graphite_CAS] = Cps
            Ts_values_dict[graphite_CAS] = Ts_Cp

    assert  298.15 in p.rawdata['T'].tolist()


# Format a file and dump it
import natsort

CASs = list(natsort.natsorted(list(set(metadata_table.index.values.tolist()))))
keys = ['CAS', 'Chemical', 'formula', 'Hfl','Gfl', 'S0l','Cpl', 'Hfg', 'Gfg', 'S0g', 'Cpg']
lines = ['\t'.join(keys) + '\n']
search_dicts = [names_dict, formulas_dict, Hfls_dict, Gfls_dict, S0ls_dict, Cpls_dict, Hfgs_dict, Gfgs_dict, S0gs_dict, Cpgs_dict]

for CAS in CASs:
    has_any = any(CAS in d for d in search_dicts[2:])
    if not has_any:
        continue
    parts = [CAS] + [str(d.get(CAS, '')) for d in search_dicts]
    line = '\t'.join(parts) + '\n'
    lines.append(line)

    # print(line)

dump_csv_path = os.path.join(dirname, '..', '..', '..', 'chemicals', 'Reactions', 'JANAF_1998.tsv')
print(dump_csv_path)
f = open(dump_csv_path, 'w')
f.writelines(lines)
f.close()

# print('Liquids found')
# print('-'*100)
# for CAS in CASs:
#     if CAS in Cpls_dict:
#         print(formulas_dict[CAS], names_dict[CAS], CAS)
#
# print('gases found')
# print('-'*100)
# for CAS in CASs:
#     if CAS in Cpgs_dict:
#         print(formulas_dict[CAS], names_dict[CAS], CAS)

to_dump_Cp_gas = {}
for CAS, Cps in Cpgs_values_dict.items():
    to_dump_Cp_gas[CAS] = (Tg_values_dict[CAS], Cps)

dump_Cpg_path = os.path.join(dirname, '..', '..', '..', 'chemicals', 'Heat Capacity', 'JANAF_1998_gas_Cp.json')
f = open(dump_Cpg_path, 'w')
json.dump(to_dump_Cp_gas, f, indent=2)
f.close()

to_dump_Cp_liq = {}
for CAS, Cps in Cpls_values_dict.items():
    to_dump_Cp_liq[CAS] = (Tl_values_dict[CAS], Cps)

dump_Cp_path = os.path.join(dirname, '..', '..', '..', 'chemicals', 'Heat Capacity', 'JANAF_1998_liq_Cp.json')
f = open(dump_Cp_path, 'w')
json.dump(to_dump_Cp_liq, f, indent=2)
f.close()


to_dump_Cp_solid = {}
for CAS, Cps in Cpss_values_dict.items():
    to_dump_Cp_solid[CAS] = (Ts_values_dict[CAS], Cps)

dump_Cp_path = os.path.join(dirname, '..', '..', '..', 'chemicals', 'Heat Capacity', 'JANAF_1998_solid_Cp.json')
f = open(dump_Cp_path, 'w')
json.dump(to_dump_Cp_solid, f, indent=2)
f.close()
