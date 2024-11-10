import pandas as pd
import requests
from io import StringIO
from chemicals.identifiers import check_CAS

# URL of the CSV file
url = "https://raw.githubusercontent.com/Gnpd/HSPiPy/refs/heads/main/hspipy/db.csv"

# Download the file
response = requests.get(url)
csv_content = StringIO(response.text)

# Read the CSV file
df = pd.read_csv(csv_content)

# Rename the Hansen solubility parameter columns
df.rename(columns={
    'D': 'HANSEN_DELTA_D',
    'P': 'HANSEN_DELTA_P',
    'H': 'HANSEN_DELTA_H'
}, inplace=True)

# Multiply Hansen parameters by 1000
for col in ['HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']:
    df[col] = df[col] * 1000

# Drop rows without a CAS
df = df.dropna(subset=['CAS'])

# Remove duplicates based on CAS
df = df.drop_duplicates(subset=['CAS'], keep='first')

# Filter CAS numbers
df = df[df['CAS'].apply(check_CAS)]

# Sort by CAS
df = df.sort_values(by='CAS')

# Select only the required columns
df = df[['CAS', 'HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']]

# Save to a .tsv file
df.to_csv("hspipy_hansen_solubility_parameters.tsv", sep='\t', index=False)

print("HSPiPy Hansen Solubility Parameters data processed and saved.")