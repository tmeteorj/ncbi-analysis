import os

from analysis.gene_location_analysis import GeneLocationAnalysis
from analysis.gene_similarity_match import GeneSimilarityMatch, OrderType
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import ExperimentConfig

# the file you want to run
gene_match_names = ['14ibsc.txt']
data_name = 'NC_000913.3.txt'

"""
match algorithm:
(1) text_distance  :  the minimum step to change one sequence to other sequence 
(2) direct_match   :　compare each gene one by one
(3) consistency    :  compare each gene one by one, and the longest matched sequence has the most top priority
(4) pattern        :  only find gene match some pattern
"""
weighted = [0, 0, 1, 0]
# only consider of the best gene in the range of candidate_distance
candidate_distance = 5
# output top 500 sequence
top_k = 500
# analysis 2 gene at the same time
batch_size = 2
# ignore mismatch with patience 2
patience = 2
continuous_mismatch_limit = 10
conditions = {
    'must': [{
        'offset': 0,
        'length': 4
    }, {
        'offset': -4,
        'length': 4
    }],
    'optional': [{
        'offset': 4,
        'length': 1
    }, {
        'offset': -5,
        'length': 1
    }]
}
ecocyc_file_name = 'Ecocyc_NC_000913.txt'
filter_sub_span = (45, 25)
output_promoter = False
order_type = OrderType.Increment
if __name__ == '__main__':
    ecocyc_file_path = os.path.join(ExperimentConfig.data_directory, ecocyc_file_name)
    for gene_match_name in gene_match_names:
        gene_path = os.path.join(ExperimentConfig.data_directory, gene_match_name)
        data_path = os.path.join(ExperimentConfig.rna_download_directory, data_name)
        similarity_match = GeneSimilarityMatch(gene_path, data_path, ExperimentConfig.output_directory,
                                               top_k=top_k,
                                               candidate_distance=candidate_distance,
                                               batch_size=batch_size,
                                               patience=patience,
                                               weighted=weighted,
                                               conditions=conditions,
                                               continuous_mismatch_limit=continuous_mismatch_limit,
                                               order_type=order_type)
        similarity_match.run()
        if OrderType.Decrement == order_type:
            gene_location_analysis = GeneLocationAnalysis(similarity_match.result_path, ecocyc_file_path,
                                                          ExperimentConfig.output_directory,
                                                          process_sub_data=weighted[2] > 0,
                                                          filter_sub_span=filter_sub_span,
                                                          output_promoter=output_promoter)
            gene_location_analysis.run()
