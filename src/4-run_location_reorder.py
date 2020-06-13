import os

from analysis.location_reorder import LocationReorder
from experiment_config import ExperimentConfig

index_file_name = 'test_index.txt'
location_file_name = 'dinQ_antisense_consistency_match_sub_location_result.txt'

if __name__ == '__main__':
    index_file_path = os.path.join(ExperimentConfig.data_directory, index_file_name)
    location_file_path = os.path.join(ExperimentConfig.output_directory, location_file_name)
    location_reorder = LocationReorder(location_file_path, index_file_path, ExperimentConfig.output_directory)
    location_reorder.run()
