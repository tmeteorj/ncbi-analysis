import os

from analysis.gene_location_analysis import GeneLocationAnalysis
from analysis.gene_similarity_match import GeneSimilarityMatch
from analysis.models.similarity_type import SimilarityType
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import ExperimentConfig

# the file you want to run
gene_match_names = ['sib_antisense.txt']
data_name = 'NC_000913.3.txt'
filter_gene_name = ''

weighted = {
    SimilarityType.TextEdit: 0,
    SimilarityType.Direct: 0,
    SimilarityType.Consistency: 1,
    SimilarityType.Pattern: 0,
    SimilarityType.Blat: 0
}
# only consider of the best gene in the range of candidate_distance
candidate_distance = 5
# output top 500 sequence
top_k = 1000
# ignore mismatch with patience 2
patience = 2
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

if __name__ == '__main__':
    ecocyc_file_path = os.path.join(ExperimentConfig.data_directory, ecocyc_file_name)
    for gene_match_name in gene_match_names:
        gene_path = os.path.join(ExperimentConfig.data_directory, gene_match_name)
        filter_gene_path = os.path.join(ExperimentConfig.data_directory, filter_gene_name) if filter_gene_name else None
        data_path = os.path.join(ExperimentConfig.rna_download_directory, data_name)

        similarity_match = GeneSimilarityMatch(gene_path, data_path, ExperimentConfig.output_directory,
                                               top_k=top_k,
                                               candidate_distance=candidate_distance,
                                               patience=patience,
                                               weighted=weighted,
                                               conditions=conditions)
        gene_location_analysis = GeneLocationAnalysis(similarity_match.result_path, ecocyc_file_path,
                                                      ExperimentConfig.output_directory,
                                                      process_sub_data=weighted.get(SimilarityType.Consistency, 0) > 0,
                                                      filter_gene_path=filter_gene_path)

        similarity_match.run(gene_location_analysis)

        gene_location_analysis.run()
