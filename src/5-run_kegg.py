import os

from analysis.kegg_analysis import KeggAnalysis
from experiment_config import ExperimentConfig

# input_file_name = '20200103.txt'
input_file_name = "kegg_genes.txt"
is_gene = True

if __name__ == '__main__':
    input_file_path = os.path.join(ExperimentConfig.data_directory, input_file_name)
    exp = KeggAnalysis(input_path=input_file_path,
                       download_directory=ExperimentConfig.kegg_download_directory,
                       output_directory=ExperimentConfig.output_directory,
                       is_gene=is_gene)
    exp.run()
