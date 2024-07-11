import requests
import pandas as pd
from tqdm import tqdm
import concurrent.futures

def get_smiles_from_name(chemical_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            return chemical_name, data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    return chemical_name, None

# Read the input CSV file
data = pd.read_csv(r"C:\Users\igorh\Documents\SOFIA_MQ\data\olfactionbase_odors_odorant_data.csv")

# Initialize an empty dictionary to store SMILES strings
smiles_dict = {}

# Define the number of workers for concurrent requests
num_workers = 40

# Use ThreadPoolExecutor to fetch SMILES in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
    futures = {executor.submit(get_smiles_from_name, name): name for name in data['Chemical name']}
    for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Converting to SMILES"):
        chemical_name, smiles = future.result()
        smiles_dict[chemical_name] = smiles

# Map the SMILES strings back to the original dataframe
data['SMILES'] = data['Chemical name'].map(smiles_dict)

# Save the dataframe with the SMILES column to a new CSV file
data.to_csv(r"C:\Users\igorh\Documents\SOFIA_MQ\data\olfactionbase_odors_odorant_data_with_smiles.csv", index=False)

print("Conversion complete. The data with SMILES has been saved.")