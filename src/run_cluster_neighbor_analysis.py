from analysis.cluster_match import ClusterMatcher
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import *

fna_name = '7a_utr_bacteria.fna'
rna_tag = '7a_utr'

if __name__ == '__main__':
    input_path = os.path.join(ExperimentConfig.data_directory, fna_name)
    cluster_match = ClusterMatcher(rna_tag, input_path, ExperimentConfig.output_directory)
    cluster_match.run()
    file_name = rna_tag + '_all_result.txt'
    input_path = os.path.join(ExperimentConfig.output_directory, file_name)
    neighbor_analysis = NeighborAnalysis(input_path, ExperimentConfig.rna_download_directory,
                                         ExperimentConfig.output_directory)
    neighbor_analysis.run()
