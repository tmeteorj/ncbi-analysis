from analysis.cluster_match import ClusterMatcher
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import *

fna_name = 'T2_sun_bacteria.fna'
rna_tag = 'T2_sun'

if __name__ == '__main__':
    input_path = os.path.join(data_directory, fna_name)
    cluster_match = ClusterMatcher(rna_tag, input_path, output_directory)
    cluster_match.run()
    file_name = rna_tag + '_all_result.txt'
    input_path = os.path.join(output_directory, file_name)
    neighbor_analysis = NeighborAnalysis(input_path, rna_download_directory, output_directory)
    neighbor_analysis.run()
