"""This file contains an exporter for the IARC carcinogen data.

The .csv file from https://monographs.iarc.fr/list-of-classifications should be
placed in the same directory as this script, and then the data file will be
updated when this is run.
"""

import pandas as pd
from natsort import natsorted

# For this case, do our best to store the worse ranking
IARC_ranking = {'1': 1000, '2A': 100, '2B': 5, '3': 0, '4': -10}
IARC_str_to_ints = {'1': 1, '2A': 11, '2B': 12, '3': 3, '4': 4}

# These are compounds who do not have a definitive assignment - but their descrption references another
# category which does have a definite assignment. Some of these CAS numbers are obsolete and not
# in the current document.
manual_ratings = {'1303-00-0': '1', '1464-53-5': '1', '2602-46-2': '1', '1937-37-7': '1', '6795-23-9': '2B',
                 '7664-93-9': '1', '10043-66-0': '1', '10098-97-2': '1', '13768-00-8': '1',
                 '14567-73-8': '1', '16071-86-6': '1', '17068-78-9': '1'}

CASRNs = []
descriptions = []
groups = []
volumes = []
years = []
df = pd.read_csv('List of Classifications – IARC Monographs on the Identification of Carcinogenic Hazards to Humans.csv')
for _, r in df.iterrows():
    CAS_str = r['CAS No.']
    Agent = r['Agent']
    Group = r['Group']
    Volume = r['Volume']
    Year = r['Year']
    Desecription = r['Additional information']
    if type(CAS_str) == str:
        for CAS in CAS_str.split(', '):
            if CAS in manual_ratings:
                Group = manual_ratings[CAS]
            skip = False
            if type(Group) != str:
                continue # skip nans - repeated entries
            try:
                existing_index = CASRNs.index(CAS)
            except:
                existing_index = -1
            if existing_index != -1:
                skip = True
                if IARC_ranking[groups[existing_index]] > IARC_ranking[Group]:
                    continue
                else:
                    CASRNs[existing_index] = CAS
                    descriptions[existing_index] = Agent
                    groups[existing_index] = Group
                    volumes[existing_index] = Volume
                    years[existing_index] = Year
            if not skip:
                CASRNs.append(CAS)
                descriptions.append(Agent)
                groups.append(Group)
                volumes.append(Volume)
                years.append(Year)
for i in range(len(groups)):
    groups[i] = IARC_str_to_ints[groups[i]]


df = pd.DataFrame({'description': descriptions, 'group': groups, 'volumes': volumes,
                  'year': years}, index=CASRNs)
df.fillna('', inplace=True)
df = df.reindex(index=natsorted(df.index))
df.index.name = 'CAS'

df.to_csv('../../chemicals/Safety/IARC Carcinogen Database.tsv', sep='\t')
