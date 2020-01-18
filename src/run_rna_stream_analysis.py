from analysis.gene_stream_analysis import GeneStreamAnalysis
from experiment_config import *

gene_file_name = 'NC_000913.3.txt'
rna_file_name = '77_ecocyc_result.txt'

if __name__ == '__main__':
    data_path = os.path.join(ExperimentConfig.rna_download_directory, gene_file_name)
    rna_path = os.path.join(ExperimentConfig.data_directory, rna_file_name)
    gene_stream_analysis = GeneStreamAnalysis(data_path, rna_path, ExperimentConfig.output_directory,
                                              mode='rna',
                                              limit=200)
    gene_stream_analysis.run()
