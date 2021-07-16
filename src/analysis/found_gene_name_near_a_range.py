import pandas as pd

from experiment_config import ExperimentConfig
from utils.gene_position_helper import GenePositionHelper
from utils.str_util import StrConverter


class FoundGeneNameNearARange:

    def __init__(self, helper: GenePositionHelper, output_directory=ExperimentConfig.output_directory):
        self.helper = helper
        self.expand_headers = ['related', 'hit', 'sequence']
        self.output_directory = output_directory

    def run(self, input_file_path):
        df = pd.read_csv(input_file_path, sep='\t')
        df[self.expand_headers] = df.apply(self.get_gene_name_and_sequence, axis=1, result_type='expand')
        output_file_path = StrConverter.generate_result_file_name(input_file_path, self.output_directory, 'near_gene')
        df.to_csv(output_file_path, sep='\t', index=False)

    def get_gene_name_and_sequence(self, record: pd.Series):
        left, right = record['end_point_1'], record['end_point_2']
        if left > right:
            left, right = right, left
        result = self.helper.get_nearby_gene_based_by_range(left, right, '+')
        return tuple([result.get(header, '') for header in self.expand_headers])
