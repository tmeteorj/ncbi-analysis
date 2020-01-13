from analysis.gene_extract import GeneExtract
from experiment_config import *

extract_gene_file_name = 'NC_000913.3.txt'
extract_gene_sequence = 'gene_utr_to_find.txt'
gene_extract_based = 'range'

if __name__ == '__main__':
    data_path = os.path.join(ExperimentConfig.rna_download_directory, extract_gene_file_name)
    rna_path = os.path.join(ExperimentConfig.data_directory, extract_gene_sequence)
    gene_extract = GeneExtract(data_path, rna_path, ExperimentConfig.output_directory, gene_extract_based)
    gene_extract.run()
