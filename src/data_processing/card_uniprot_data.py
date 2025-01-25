import os
import pandas as pd
from Bio import SeqIO
from io import StringIO
import gzip
import re


# Directories
BASE_DIR_CARD = "./data/raw/card-data"
BASE_DIR_UNIPROT = "./data/raw/uniprot"
SAVE_DIR = "./data/processed"


def load_card_data(base_dir):
    file_names = os.listdir(base_dir)
    meta_files = [f for f in file_names if f.endswith('.tsv')]
    
    if len(meta_files) < 1:
        raise FileNotFoundError("Expected at least 1 metadata files in the CARD directory.")
    
    card_dict = {}
    for m in meta_files:
        card_dict[m.split('.tsv')[0]] = pd.read_csv(f"{base_dir}/{m}", sep='\t')
    
    return card_dict

# tmp = load_card_data(BASE_DIR_CARD)
# print(tmp.keys())


def parse_fasta_file(fasta_path):
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append({
            "id": record.id,
            "description": record.description,
            "sequence": str(record.seq)
        })
    
    return pd.DataFrame(sequences)


def clean_string(data_str, category):
    # aro_categories_index, species
    if category.lower() == 'species':
        result = data_str.split('[')[-1]
        return result.split(']')[0]
    # uniprot, aro
    if category.lower() == 'aro':
        target = r"GN=[a-zA-Z0-9]+"
        x = re.rearch(target, data_str)
        return x.group()
    

def clean_card_data(seq_df, aro_index):
    seq_df['ARO Accession'] = seq_df['id'].apply(lambda x: x.split("|")[-2])
    seq_df.set_index('ARO Accession', inplace=True)
    aro_index.set_index('ARO Accession', inplace=True)
    
    merged = seq_df.join(aro_index, on='ARO Accession')
    merged['species'] = merged['description'].apply(lambda x: clean_string(x, 'species'))
    merged.reset_index(inplace=True)
    # print(f"Before merging: {seq_df.shape}, After merging: {merged.shape}") # (6048, 3), (6050, 16)
    
    # Check duplicated data
    merged['duplicated'] = merged['ARO Accession'].duplicated()
    # print(merged[merged['duplicated']==True])
    merged.drop('duplicated', axis=1, inplace=True)
    
    # Remove duplicates
    duplicates_dropped = merged.drop_duplicates(subset='ARO Accession', keep='first')
    # print(f"Before removing duplcates: {merged.shape}, After removing duplicates: {duplicates_dropped.shape}") # (6050, 16), (6048, 16)
    
    return duplicates_dropped


def parse_uniprot_fast(file_path):
    with gzip.open(file_path, "rt") as f:
        fasta_data =f.read()
        
    sequences = []
    for record in SeqIO.parse(StringIO(fasta_data), "fasta"):
        sequences.append({
            "id": record.id,
            "description": record.description,
            "sequence": str(record.seq)
        })
    
    return pd.DataFrame(sequences) 


def resistance_check(description, keywords):
    return any(keyword.lower() in description.lower() for keyword in keywords)


def filter_non_resistance_data(non_res_df, card_data, keywords):
    # Filter with resistance keywords
    non_res_df['suspicious'] = non_res_df['description'].apply(lambda x: resistance_check(x, keywords))
    # Filter with ARO Name from CARD data
    card_keywords = [aro.replace("-", "_") for aro in list(card_data['ARO Name'])]
    non_res_df['suspicious'] = non_res_df['description'].apply(lambda x: resistance_check(x, card_keywords))
    filtered = non_res_df[non_res_df['suspicious'] == False]
    print(f"Shape of filtered Uniprot data: {filtered.shape}")
    
    return filtered


if __name__ == "__main__":
    # Load CARD meta-data
    card_dict = load_card_data(BASE_DIR_CARD)    
    
    # Parse CARD FASTA file
    seq_df = parse_fasta_file(f"{BASE_DIR_CARD}/protein_fasta_protein_homolog_model.fasta")
    merged_card_data = clean_card_data(seq_df, card_dict['aro_index'])
    merged_card_data.to_csv(f"{SAVE_DIR}/merged_aro.csv")
    
    # Parse UniProt FASTA data (manually downloaded from UniProt)
    uniprot_file_path = f"{BASE_DIR_UNIPROT}/uniprotkb_bacteria_NOT_resistance_beta_2025_01_16.fasta.gz"
    uniprot_sequences = parse_uniprot_fast(uniprot_file_path)
    
    # Filter non-resistance data from uniprot sequences
    resistance_keywords = ["resistance", "beta-lactamase", "aminoglycoside", "macrolide", "tetracycline"]
    filtered_uniprot = filter_non_resistance_data(uniprot_sequences, merged_card_data, resistance_keywords)
    filtered_uniprot.to_csv(f"{SAVE_DIR}/filtered_uniprot.csv")
    
