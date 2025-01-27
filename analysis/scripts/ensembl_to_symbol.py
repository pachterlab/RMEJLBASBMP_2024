import argparse
import csv
from tqdm import tqdm
from gget.gget_info import info

def ensembl_to_gene_name_for_single_gene(gene_id):
    """
    Function to fetch gene name from a single Ensembl IDs string using gget info.
    """
    # Remove version number if passed
    gene_id = gene_id.split(".")[0]

    try:
        info_df = info(gene_id, pdb=False, ncbi=False, uniprot=False, verbose=False)

        # Check if Ensembl ID was found
        if isinstance(info_df, type(None)):
            return None

        gene_symbol = info_df.loc[gene_id]["ensembl_gene_name"]
    
        # If more than one gene symbol was returned, use first entry
        if isinstance(gene_symbol, list):
            gene_symbol = str(gene_symbol[0])
        else:
            gene_symbol = str(gene_symbol)
    except Exception as e:
        gene_symbol = None

    return gene_symbol

def main():
    # define unique_genes_py in argparse as a positional argument
    parser = argparse.ArgumentParser(description="Convert Ensembl IDs to gene symbols")
    parser.add_argument("-i", "--input_gene_list_file_path", type=str, help="List of unique Ensembl IDs")
    parser.add_argument("-o", "--output_csv_path", type=str, help="Path to save the gene dictionary as a CSV")
    args = parser.parse_args()

        # Load the text file as a list
    with open(args.input_gene_list_file_path, "r") as file:
        unique_genes_py = [line.strip() for line in file]

    # Initialize an empty list for the results
    unique_gene_symbols_py = []

    for gene_id in tqdm(unique_genes_py, desc="Processing genes"):
        unique_gene_symbols_py.append(ensembl_to_gene_name_for_single_gene(gene_id))

    gene_dict = dict(zip(unique_genes_py, unique_gene_symbols_py))

    # Save the dictionary as a CSV
    with open(args.output_csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Ensembl_ID", "Gene_Symbol"])  # Write the header
        for ensembl_id, gene_symbol in gene_dict.items():
            writer.writerow([ensembl_id, gene_symbol])

if __name__ == "__main__":
    main()