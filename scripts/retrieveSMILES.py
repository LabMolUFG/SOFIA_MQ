import pubchempy as pcp
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

def get_smiles_from_name(chemical_name, retries=5):
    for i in range(retries):
        try:
            compound = pcp.get_compounds(chemical_name, 'name')
            if compound:
                return compound[0].canonical_smiles
        except pcp.PubChemHTTPError as e:
            if 'PUGREST.ServerBusy' in str(e):
                wait_time = 2 ** i  # Exponential backoff
                time.sleep(wait_time)
            else:
                print(f"Error: {e}")
                break
    return None

# Read the input CSV file
data = pd.read_csv(r"C:\Users\igorh\Documents\SOFIA_MQ\data\olfactionbase_odors_odorant_data.csv")

# Function to process each chemical name
def process_chemical_name(chemical_name):
    return chemical_name, get_smiles_from_name(chemical_name)

# Initialize an empty list with None values to store SMILES strings
smiles_list = [None] * len(data)

# Use ThreadPoolExecutor for multithreading
with ThreadPoolExecutor(max_workers=10) as executor:
    # Submit tasks to the executor and keep track of the futures
    futures = {executor.submit(process_chemical_name, name): idx for idx, name in enumerate(data['Chemical name'])}
    
    for future in tqdm(as_completed(futures), total=len(futures), desc="Converting to SMILES"):
        idx = futures[future]
        try:
            chemical_name, smiles = future.result()
            smiles_list[idx] = smiles
        except Exception as e:
            print(f"Error processing {data['Chemical name'][idx]}: {e}")

# Add the SMILES list to the dataframe
data['SMILES'] = smiles_list

# Save the dataframe with the SMILES column to a new CSV file
data.to_csv(r"C:\Users\igorh\Documents\SOFIA_MQ\data\olfactionbase_odors_odorant_data_with_smiles.csv", index=False)

print("Conversion complete. The data with SMILES has been saved.")
