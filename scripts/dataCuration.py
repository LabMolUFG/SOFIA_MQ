import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Definition of odor column
odorcolumn = "Sub-odor"

def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Standardize the molecule
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)
    return Chem.MolToSmiles(mol)

def process_smiles(input_csv, removed_csv):
    df = pd.read_csv(input_csv)
    if 'SMILES' not in df.columns or odorcolumn not in df.columns:
        raise ValueError(f"The input CSV file must contain 'SMILES' and '{odorcolumn}' columns")
    
    # Remove rows with None values
    df.dropna(inplace=True)

    # Separate SMILES that represent mixtures
    mixed_df = df[df['SMILES'].str.contains(r'\.')]
    mixed_df.to_csv(removed_csv, index=False)

    # Remove SMILES that represent mixtures from the main DataFrame
    df = df[~df['SMILES'].str.contains(r'\.')]

    df['standardized_smiles'] = df['SMILES'].apply(standardize_smiles)
    
    # Group by standardized SMILES and aggregate odors
    grouped = df.groupby('standardized_smiles')[odorcolumn].apply(lambda x: '/'.join(x)).reset_index()
    
    return grouped

# Function to handle splitting columns with "/" and strip spaces
def split_columns(df):
    # Iterate over columns starting from the second column (index 1)
    for col in df.columns[1:]:
        # Split values in each column by "/" and strip spaces
        df[col] = df[col].str.split('/').apply(lambda x: [i.strip() for i in x] if isinstance(x, list) else x)
    return df

# Function to filter out corrupted labels
def filter_corrupted_labels(odors):
    return [odor for odor in odors if len(odor) > 1 and not odor.endswith('&')]

# Function to remove duplicates and redistribute odors into separate columns
def remove_and_redistribute_duplicates(df):
    new_rows = []
    
    for _, row in df.iterrows():
        unique_odors = set()
        for col in row[1:]:
            if isinstance(col, list):
                # Filter out corrupted labels
                col = filter_corrupted_labels(col)
                for odor in col:
                    if odor not in unique_odors:
                        unique_odors.add(odor)
        
        new_row = [row[0]] + list(unique_odors)
        new_rows.append(new_row)
    
    # Determine the maximum number of columns needed
    max_columns = max(len(row) for row in new_rows)
    column_names = ['standardized_smiles'] + [f'odor_{i+1}' for i in range(1, max_columns)]
    
    # Create a new DataFrame
    new_df = pd.DataFrame(new_rows, columns=column_names).fillna('')
    
    return new_df

if __name__ == "__main__":
    input_csv = r'C:\Users\igorh\Documents\SOFIA_MQ\data\olfactionbase_odors_odorant_data_with_smiles.csv'  # Path to the input CSV file
    removed_csv = r'C:\Users\igorh\Documents\SOFIA_MQ\data\removed_mixtures.csv'  # Path to the output CSV file for removed mixtures
    grouped = process_smiles(input_csv, removed_csv)

    # Split labels separated by "/" and strip spaces
    splitted = split_columns(grouped)
    
    # Apply the function to remove duplicates and redistribute odors
    cleaned_df = remove_and_redistribute_duplicates(splitted)
    
    # Save the cleaned DataFrame to a new CSV file
    cleaned_df.to_csv(r'C:\Users\igorh\Documents\SOFIA_MQ\data\curated_SubOdor.csv', index=False)
