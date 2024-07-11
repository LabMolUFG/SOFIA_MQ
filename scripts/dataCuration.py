import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolStandardize

# Definição da coluna de odor
odorcolumn = "Sub-odor"

def standardize_smiles(smiles):
    normalizer = MolStandardize.normalize.Normalizer()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = normalizer.normalize(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def process_smiles(input_csv, output_csv, removed_csv):
    df = pd.read_csv(input_csv, sep=";")
    if 'SMILES' not in df.columns or odorcolumn not in df.columns:
        raise ValueError(f"The input CSV file must contain 'SMILES' and '{odorcolumn}' columns")

    # Separar SMILES que representam misturas
    mixed_df = df[df['SMILES'].str.contains(r'\.')]
    mixed_df.to_csv(removed_csv, index=False)

    # Remover SMILES que representam misturas do DataFrame principal
    df = df[~df['SMILES'].str.contains(r'\.')]

    df['standardized_smiles'] = df['SMILES'].apply(standardize_smiles)
    
    # Agrupar por SMILES padronizadas e agregar odores
    grouped = df.groupby('standardized_smiles')[odorcolumn].apply(lambda x: x.tolist()).reset_index()
    
    # Criar um DataFrame com colunas separadas para cada odor
    max_odors = grouped[odorcolumn].apply(len).max()
    for i in range(max_odors):
        grouped[f'odor_{i+1}'] = grouped[odorcolumn].apply(lambda x: x[i] if i < len(x) else None)
    
    grouped.drop(columns=[odorcolumn], inplace=True)
    
    grouped.to_csv(output_csv, index=False)

if __name__ == "__main__":
    input_csv = 'teste.csv'  # Caminho para o arquivo CSV de entrada
    output_csv = 'standardized_smiles_with_odors.csv'  # Caminho para o arquivo CSV de saída
    removed_csv = 'removed_mixtures.csv'  # Caminho para o arquivo CSV de saída das misturas removidas

    process_smiles(input_csv, output_csv, removed_csv)
