import os

from analysis.gene_location_analysis import format_data_to_tsv
from experiment_config import ExperimentConfig
from utils.ecocyc_data_loader import EcocycDataLoader

ecocyc_file_name = 'Ecocyc_NC_000913_2.0.txt'
input_file_names = ["mgtS_antisense_consistency_match_sub_location_result.txt"]

if __name__ == '__main__':
    ecocyc_file_path = os.path.join(ExperimentConfig.data_directory, ecocyc_file_name)
    ecocyc_data_loader = EcocycDataLoader(ecocyc_file_path)
    ecocyc_data_loader.build_database()
    for file_name in input_file_names:
        file_path = os.path.join(ExperimentConfig.output_directory, file_name)
        format_data_to_tsv(file_path, file_path.replace('.txt', '_format.tsv'), ecocyc_data_loader)
