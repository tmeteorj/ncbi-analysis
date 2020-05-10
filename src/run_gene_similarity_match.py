import os

from analysis.gene_similarity_match import GeneSimilarityMatch
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import ExperimentConfig

gene_match_names = ['dinQ.txt']
data_name = 'NC_000913.3.txt'
match_algorithm = 'text_distance'
# match_algorithm = 'char_match'
candidate_distance = 5
top_k = 30
batch_size = 2

if __name__ == '__main__':
    for gene_match_name in gene_match_names:
        gene_path = os.path.join(ExperimentConfig.data_directory, gene_match_name)
        data_path = os.path.join(ExperimentConfig.rna_download_directory, data_name)
        similarity_match = GeneSimilarityMatch(gene_path, data_path, ExperimentConfig.output_directory,
                                               top_k=top_k,
                                               match_algorithm=match_algorithm,
                                               candidate_distance=candidate_distance,
                                               batch_size=batch_size)
        similarity_match.run()

        neighbor_analysis = NeighborAnalysis(similarity_match.result_path, ExperimentConfig.rna_download_directory,
                                             ExperimentConfig.output_directory)
        neighbor_analysis.run()