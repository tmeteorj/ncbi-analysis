import json
import re

import pandas as pd
from experiment_config import ExperimentConfig
from utils.coli_database import ColiDatabase
from utils.str_util import StrConverter


class GentamycinAnalysis:
    __header_gene__ = 'Gene'
    __expand_headers__ = ['left_gene', 'right_gene', 'hit_gene', 'sequence']

    def __init__(self, coli_database_path, output_directory=None):
        self.coli_database = ColiDatabase(coli_database_path)
        self.output_directory = output_directory if output_directory else ExperimentConfig.output_directory

    def run(self, gene_list_path: str):
        output_path = StrConverter.generate_result_file_name(gene_list_path, self.output_directory, 'gentamycin')
        gene_df = pd.read_csv(gene_list_path, sep='\t')
        gene_df[self.__expand_headers__] = gene_df.apply(lambda record: self.expand_one_record(record), axis=1,
                                                         result_type='expand')
        gene_df.to_csv(output_path, sep='\t', index=None, header=True)

    def expand_one_record(self, record: pd.Series):
        """

        :param record:
        :return: left gene, right gene, hit gene, sequence
        """
        if not record[self.__header_gene__].startswith('DR'):
            result = tuple(['' for _ in self.__expand_headers__])
            return result
        else:
            result = {}
            location_column = record['Locus'].strip()
            matched = re.findall(r'(.+):(\d+)-(\d+)\((.)\)', location_column)
            assert matched and len(matched) == 1, location_column
            matched = matched[0]
            left, right, direction = int(matched[1]), int(matched[2]), matched[3]
            left_ge_id = self.coli_database.find_first_greater_equal(left)
            left_lt_id = left_ge_id - 1
            right_ge_id = self.coli_database.find_first_greater_equal(right)
            right_lt_id = right_ge_id - 1
            if left_ge_id == right_lt_id:
                result['hit_gene'] = self.coli_database.segments[left_ge_id].gene
                result['sequence'] = self.coli_database.segments[left_ge_id].sequence
            elif left_ge_id < right_lt_id:
                assert left_ge_id + 1 == right_lt_id
                result['left_gene'] = self.coli_database.segments[left_ge_id].gene
                result['right_gene'] = self.coli_database.segments[right_lt_id].gene
                result['sequence'] = 'left:' + self.coli_database.segments[left_ge_id].sequence + \
                                     ',right：' + self.coli_database.segments[right_lt_id].sequence
            else:
                assert left_ge_id - 1 == right_lt_id
                result['left_gene'] = self.coli_database.segments[right_lt_id].gene
                result['right_gene'] = self.coli_database.segments[left_ge_id].gene
                result['sequence'] = 'left:' + self.coli_database.segments[right_lt_id].sequence + \
                                     ',right:' + self.coli_database.segments[left_ge_id].sequence
            result = tuple([result.get(header, '') for header in self.__expand_headers__])
            return result
