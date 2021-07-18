import pandas as pd

from analysis.gentamycin import GentamycinAnalysis
from experiment_config import ExperimentConfig
from utils.gene_position_helper import GenePositionHelper
from utils.str_util import StrConverter


class FoundGeneNameNearARange:
    __keep_headers__ = ['name', 'file', 'length', 'start', 'end', 'compare_length', 'type', 'gene', 'gene_left',
                        'gene_right', 'sequence']

    def __init__(self, helper: GenePositionHelper, output_directory=ExperimentConfig.output_directory):
        self.helper = helper
        self.expand_headers = ['related', 'hit', 'sequence']
        self.output_directory = output_directory

    def run(self, input_file_path):
        df = pd.read_csv(input_file_path, sep='\t')
        df[['found', 'length', 'compare_length']] = df.apply(self.get_gene_name_and_sequence, axis=1,
                                                             result_type='expand')
        explode_list = []
        for _, record in df.iterrows():
            explode_result = {column: record.get(column, '') for column in self.__keep_headers__}
            for gene_record in record['found']:
                new_result = {k: v for k, v in explode_result.items()}
                new_result.update({
                    k: v for k, v in gene_record.items()
                })
                explode_list.append(new_result)
        df = pd.DataFrame(explode_list)
        output_file_path = StrConverter.generate_result_file_name(input_file_path, self.output_directory, 'near_gene')
        df[self.__keep_headers__].to_csv(output_file_path, sep='\t', index=False)

    def get_gene_name_and_sequence(self, record: pd.Series):
        start, end = record['start'], record['end']
        if start > end:
            left, right = end, start
            direction = '-'
        else:
            left, right = start, end
            direction = '+'
        results = list(self.helper.get_nearby_gene_based_by_range(left, right, direction))
        compare_length = right - left + 1
        left, right, direction = GentamycinAnalysis.get_position(record['locus'])
        return results, right - left + 1, compare_length
