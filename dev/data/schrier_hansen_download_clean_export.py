import pandas as pd
import requests
import re
from io import StringIO
from chemicals.identifiers import CAS_from_any, check_CAS

# URL of the TSV file
url = "https://www.wolframcloud.com/obj/fcd55d96-92eb-43bb-9032-0eb2da9ddc4f"

# Download the file
response = requests.get(url)
response.encoding = 'utf-8'  # Specify UTF-8 encoding
tsv_content = response.text
# Read the TSV file
df = pd.read_csv(StringIO(tsv_content), sep='\t')

# Function to clean chemical entity and extract name
def clean_chemical(entity):
    match = re.search(r'Entity\["Chemical",\s*"([^"]+)"\]', entity)
    return match.group(1) if match else entity

# Function to extract numerical value from "Quantity[...]"
def extract_quantity(quantity_str):
    match = re.search(r'Quantity\[(\d+(\.\d+)?),', quantity_str)
    return float(match.group(1)) if match else None

# Function to find CAS number by InChIKey or compound name
def find_cas(inchi_key, compound_name):
    try:
        cas_number = CAS_from_any('InChIKey=' + inchi_key)
        return cas_number
    except ValueError:
        try:
            cas_number = CAS_from_any(compound_name)
            return cas_number
        except ValueError:
            print(f'Failed to find CAS for: InChIKey={inchi_key}, Compound={compound_name}')
            return None

# Clean up the 'Solvent' column to extract chemical names
df['Solvent'] = df['Solvent'].apply(clean_chemical)

# Extract numerical values from columns with "Quantity[...]"
df['HANSEN_DELTA_D'] = df['δd'].apply(extract_quantity)
df['HANSEN_DELTA_P'] = df['δp'].apply(extract_quantity)
df['HANSEN_DELTA_H'] = df['δh'].apply(extract_quantity)

# Multiply Hansen parameters by 1000
for col in ['HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']:
    df[col] = df[col] * 1000

# Drop the original δd, δp, δh columns
df.drop(columns=['δd', 'δp', 'δh', 'Volume', 'δt'], inplace=True)

# Drop rows without Hansen parameters
df = df.dropna(subset=['HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H'])

# Find the CAS number using InChIKey and compound name
df['CAS'] = df.apply(lambda row: find_cas(row['InChIKey'], row['Solvent']), axis=1)

# Drop rows without CAS numbers
df = df.dropna(subset=['CAS'])

# Remove duplicates based on CAS
df = df.drop_duplicates(subset=['CAS'], keep='first')

# Filter CAS numbers
df = df[df['CAS'].apply(check_CAS)]

# Sort by CAS for consistency
df = df.sort_values(by='CAS')

# Select relevant columns for the final output
df = df[['CAS', 'HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']]

# Save to a .tsv file
df.to_csv("schrier_hansen_solubility_parameters.tsv", sep='\t', index=False)

print("Joshua Schrier Hansen Solubility Parameters data processed and saved.")