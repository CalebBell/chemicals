import pandas as pd
import requests
from io import BytesIO
from chemicals.identifiers import CAS_from_any, check_CAS

# URL of the Excel file
url = "https://www.hansen-solubility.com/contents/HSP_Calculations_Ruben-Manuel.xlsx"

# Download the file
response = requests.get(url)
excel_file = BytesIO(response.content)

# Read the Excel file, specifically the 'Chemicals' sheet
df = pd.read_excel(excel_file, sheet_name='Chemicals')

# Function to validate and find CAS numbers
def validate_and_find_cas(cas, compound_name):
    if check_CAS(cas):
        return cas
    else:
        try:
            return CAS_from_any(compound_name)
        except:
            print(f'Failed to find CAS for: Compound={compound_name}')
            return None

# Apply the function to validate or find CAS numbers
df['CAS'] = df.apply(lambda row: validate_and_find_cas(row['CAS'], row['Chemical']), axis=1)

# Drop rows without valid CAS numbers
df = df.dropna(subset=['CAS'])

# Rename the Hansen solubility parameter columns
df.rename(columns={
    'dD': 'HANSEN_DELTA_D',
    'dP': 'HANSEN_DELTA_P',
    'dH': 'HANSEN_DELTA_H'
}, inplace=True)

# Multiply Hansen parameters by 1000
for col in ['HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']:
    df[col] = df[col] * 1000

# Remove duplicates based on CAS
df = df.drop_duplicates(subset=['CAS'], keep='first')

# Filter CAS numbers
df = df[df['CAS'].apply(check_CAS)]

# Sort by CAS
df = df.sort_values(by='CAS')

# Select relevant columns for the final output
df = df[['CAS', 'HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']]

# Save to a .tsv file
df.to_csv("ruben_manuel_hansen_solubility_parameters.tsv", sep='\t', index=False)

print("HSP Calculations Ruben-Manuel data processed and saved.")