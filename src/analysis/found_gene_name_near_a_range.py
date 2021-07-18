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
        expode = []
        for _, record in df.iterrows():
            expode_result = {column: record[column] for column in record.keys()}
            if isinstance(record['related'], list):
                if len(record['related']) > 1:
                    for index, (gene, sequence) in enumerate(zip(record['related'], record['sequence']['related'])):
                        new_result = {k: v for k, v in expode_result.items()}
                        new_result['related'] = str(index + 1) + '-' + gene
                        new_result['sequence'] = sequence
                        expode.append(new_result)
                else:
                    new_result = {k: v for k, v in expode_result.items()}
                    new_result['related'] = record['related'][0]
                    new_result['sequence'] = record['sequence']['related'][0]
                    expode.append(new_result)
            else:
                expode_result['sequence'] = record['sequence']['hit']
                expode.append(expode_result)
        df = pd.DataFrame(expode)
        output_file_path = StrConverter.generate_result_file_name(input_file_path, self.output_directory, 'near_gene')
        df.to_csv(output_file_path, sep='\t', index=False)

    def get_gene_name_and_sequence(self, record: pd.Series):
        start, end = record['start'], record['end']
        if start > end:
            left, right = end, start
            direction = '-'
        else:
            left, right = start, end
            direction = '+'
        result = self.helper.get_nearby_gene_based_by_range(left, right, direction)
        return tuple([result.get(header, '') for header in self.expand_headers])
