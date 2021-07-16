import os

from analysis.gentamycin import GentamycinAnalysis
from experiment_config import ExperimentConfig

gene_list_files = ['gentamycin\\down.txt', 'gentamycin\\up.txt']
atcc_database_path = os.path.join(ExperimentConfig.rna_download_directory, 'CP009072.1.txt')
ncbi_database_path = os.path.join(ExperimentConfig.rna_download_directory, 'sequence.gb')

if __name__ == '__main__':
    analysis = GentamycinAnalysis(ncbi_database_path=ncbi_database_path)
    for gene_list_file in gene_list_files:
        analysis.run(os.path.join(ExperimentConfig.data_directory, gene_list_file))
