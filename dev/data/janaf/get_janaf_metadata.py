import os

import pandas as pd

# download all the data
#os.system('wget -m https://janaf.nist.gov/ -e robots=off')

dirname = os.path.dirname(os.path.abspath(__file__))


base_folder = os.path.join(dirname, 'janaf.nist.gov', 'tables')
index_files = os.listdir(base_folder)
index_files = [i for i in index_files if i.endswith('-index.html')]
index_files = [os.path.join(base_folder, n) for n in index_files]
all_tables = []
for file in index_files:
    table = pd.read_html(file)[0]
#     print(table.shape)
    all_tables.append(table)

metadata_table = pd.concat(all_tables, axis=0)
metadata_table = metadata_table.drop_duplicates()

metadata_table = metadata_table.set_axis(metadata_table.iloc[0], axis=1)

# drop the first row
metadata_table = metadata_table.iloc[1:]
# drop two others
metadata_table.drop(inplace=True, columns=['JANAF Table', 'Links'])
# set the index
metadata_table.set_index('CAS Number', inplace=True)


metadata_table.to_csv(os.path.join(dirname, 'janaf_tables_index.csv'))


