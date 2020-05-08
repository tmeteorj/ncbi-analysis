from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import *

file_name = 'gene_14_match_result.txt'

if __name__ == '__main__':
    input_path = os.path.join(ExperimentConfig.output_directory, file_name)
    neighbor_analysis = NeighborAnalysis(input_path, ExperimentConfig.rna_download_directory,
                                         ExperimentConfig.output_directory)
    neighbor_analysis.run()
