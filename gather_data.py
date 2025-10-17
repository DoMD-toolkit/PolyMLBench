import os

import numpy as np
import pandas as pd

# rdkit descriptors
# Specify the directory you want to walk through
directory = 'data'

data_dic = {}
all_prop = []
# Walk through the directory
for dirpath, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        # Construct the full file path
        file_path = os.path.join(dirpath, filename)
        if file_path.endswith('.csv'):
            property_name = filename.replace('.csv', '')
            all_prop.append(property_name)

all_rows = 0
for dirpath, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        # Construct the full file path
        file_path = os.path.join(dirpath, filename)
        if file_path.endswith('.csv'):
            property_name = filename.replace('.csv', '')
            data = pd.read_csv(file_path)
            for index, row in data.iterrows():
                if data_dic.get(row['smiles']) is None:
                    data_dic[row['smiles']] = {}
                    for property_name in all_prop:
                        data_dic[row['smiles']][property_name] = np.nan
                data_dic[row['smiles']][property_name] = row['value']

df = pd.DataFrame.from_dict(data_dic, orient='index')

df.reset_index(inplace=True)
df.rename(columns={'index': 'smiles'}, inplace=True)
df.to_csv('data.csv', index=False)
