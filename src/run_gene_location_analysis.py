import os

from analysis.gene_location_analysis import GeneLocationAnalysis
from experiment_config import ExperimentConfig

location_analysis_names = ['dsrA_antisense_consistency_match_result.txt']
ecocyc_file_name = 'Ecocyc_NC_000913.txt'

if __name__ == '__main__':
    ecocyc_file_path = os.path.join(ExperimentConfig.data_directory, ecocyc_file_name)
    for input_file_name in location_analysis_names:
        input_file_path = os.path.join(ExperimentConfig.data_directory, input_file_name)
        gene_location_analysis = GeneLocationAnalysis(input_file_path, ecocyc_file_path,
                                                      ExperimentConfig.output_directory)
        gene_location_analysis.run()
