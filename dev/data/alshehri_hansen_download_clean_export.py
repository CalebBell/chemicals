import pandas as pd
import requests
from io import BytesIO
from chemicals.identifiers import CAS_from_any, check_CAS
from rdkit import Chem

def get_cas(smiles, compound_name):
    try:
        cas_number = CAS_from_any(smiles)
        return cas_number
    except ValueError:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            inchi_key = Chem.MolToInchiKey(mol)
            try:
                return CAS_from_any('InChIKey=' + inchi_key)
            except ValueError:
                try:
                    return CAS_from_any(compound_name)
                except ValueError:
                    print(f'Failed to get CAS for: SMILES={smiles}, InChIKey={inchi_key}, Compound={compound_name}')
                    return None
        else:
            try:
                return CAS_from_any(compound_name)
            except ValueError:
                print(f'Failed to get CAS for: Compound={compound_name}')
                return None

# URL of the Excel file
url = "https://github.com/PEESEgroup/Pure-Component-Property-Estimation/raw/main/2%20Full%20Excel%20Sheet/Full%20Results.xlsx"

# Download the file
response = requests.get(url)
excel_file = BytesIO(response.content)

# Read the Excel sheets
df_D = pd.read_excel(excel_file, sheet_name='Delta_D')
df_P = pd.read_excel(excel_file, sheet_name='Delta_P')
df_H = pd.read_excel(excel_file, sheet_name='Delta_H')

# Rename columns
df_D.rename(columns={'Experimental': 'Delta_D'}, inplace=True)
df_P.rename(columns={'Experimental': 'Delta_P'}, inplace=True)
df_H.rename(columns={'Experimental': 'Delta_H'}, inplace=True)

# Merge DataFrames
merged_df = df_D.merge(df_P, on='SMILES').merge(df_H, on='SMILES')

# Drop rows without all three values
final_df = merged_df.dropna(subset=['Delta_D', 'Delta_P', 'Delta_H'])

# Add CAS column
final_df['CAS'] = final_df.apply(lambda row: get_cas(row['SMILES'], row['Compound']), axis=1)

# Drop rows without CAS
final_df = final_df.dropna(subset=['CAS'])

# Remove duplicates based on CAS
final_df = final_df.drop_duplicates(subset=['CAS'], keep='first')

# Sort by CAS
final_df = final_df.sort_values(by='CAS')

# Rename columns and multiply values by 1000
final_df.rename(columns={
    'Delta_D': 'HANSEN_DELTA_D',
    'Delta_P': 'HANSEN_DELTA_P',
    'Delta_H': 'HANSEN_DELTA_H'
}, inplace=True)

for col in ['HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']:
    final_df[col] = final_df[col] * 1000

# Select only the required columns
final_df = final_df[['CAS', 'HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']]

# Filter CAS numbers
final_df = final_df[final_df['CAS'].apply(check_CAS)]

# Save to a .tsv file
final_df.to_csv("alshehri_hansen_solubility_parameters.tsv", sep='\t', index=False)

print("Alshehri Hansen Solubility Parameters data processed and saved.")