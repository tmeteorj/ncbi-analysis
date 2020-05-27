import os

from analysis.gene_location_analysis import GeneLocationAnalysis
from analysis.gene_similarity_match import GeneSimilarityMatch
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import ExperimentConfig

# the file you want to run
gene_match_names = ['dsrA_antisense.txt']
data_name = 'NC_000913.3.txt'

"""
match algorithm:
(1) text_distance  :  the minimum step to change one sequence to other sequence 
(2) char_match     :　compare each gene one by one
(3) consistency    :  compare each gene one by one, and the longest matched sequence has the most top priority
"""
# match_algorithm = 'text_distance'
# match_algorithm = 'char_match'
match_algorithm = 'consistency'

# only consider of the best gene in the range of candidate_distance
candidate_distance = 5
# output top 500 sequence
top_k = 500
# analysis 2 gene at the same time
batch_size = 2
# ignore gene sequence whose similarity is less than min_similarity
min_similarity = 0.01
# ignore mismatch with patience 2
patience = 2

ecocyc_file_name = 'Ecocyc_NC_000913.txt'

if __name__ == '__main__':
    ecocyc_file_path = os.path.join(ExperimentConfig.data_directory, ecocyc_file_name)
    for gene_match_name in gene_match_names:
        gene_path = os.path.join(ExperimentConfig.data_directory, gene_match_name)
        data_path = os.path.join(ExperimentConfig.rna_download_directory, data_name)
        similarity_match = GeneSimilarityMatch(gene_path, data_path, ExperimentConfig.output_directory,
                                               top_k=top_k,
                                               match_algorithm=match_algorithm,
                                               candidate_distance=candidate_distance,
                                               batch_size=batch_size,
                                               min_similarity=min_similarity,
                                               patience=patience)
        similarity_match.run()

        gene_location_analysis = GeneLocationAnalysis(similarity_match.result_path, ecocyc_file_path,
                                                      ExperimentConfig.output_directory)
        gene_location_analysis.run()