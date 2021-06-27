import os

from analysis.gene_location_analysis import GeneLocationAnalysis
from experiment_config import ExperimentConfig

location_analysis_names = ['14ibsc_match_result.txt']
ecocyc_file_name = 'Ecocyc_NC_000913.txt'
filter_sub_span = (45, 25)
output_promoter = False
if __name__ == '__main__':
    ecocyc_file_path = os.path.join(ExperimentConfig.data_directory, ecocyc_file_name)
    for input_file_name in location_analysis_names:
        input_file_path = os.path.join(ExperimentConfig.output_directory, input_file_name)
        gene_location_analysis = GeneLocationAnalysis(input_file_path,
                                                      ecocyc_file_path,
                                                      ExperimentConfig.output_directory,
                                                      process_sub_data=True,
                                                      filter_sub_span=filter_sub_span)
        gene_location_analysis.run()
