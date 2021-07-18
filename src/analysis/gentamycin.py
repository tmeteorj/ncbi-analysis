import json
import re

import pandas as pd
from experiment_config import ExperimentConfig
from utils.atcc_database import ATCCDatabase
from utils.gene_position_helper import GenePositionHelper
from utils.ncbi_database import NCBIDatabase
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter


class GentamycinAnalysis:
    __header_gene__ = 'Gene'
    __header_sequence__ = 'sequence'
    __atcc_expand_headers__ = ['related', 'hit', 'sequence']
    __ncbi_expand_headers__ = ['sequence']

    def __init__(self, ncbi_database_path=None, atcc_database_path=None, output_directory=None):
        self.atcc_database = ATCCDatabase(atcc_database_path) if atcc_database_path else None
        self.ncbi_database = NCBIDatabase(ncbi_database_path) if ncbi_database_path else None
        if self.ncbi_database:
            self.gene_position_helper = GenePositionHelper(self.ncbi_database)
        else:
            self.gene_position_helper = GenePositionHelper(self.atcc_database)
        assert self.atcc_database is not None or self.ncbi_database is not None
        self.output_directory = output_directory if output_directory else ExperimentConfig.output_directory
        self.expand_headers = self.__atcc_expand_headers__ if self.atcc_database else self.__ncbi_expand_headers__
        self.expand_one_record_func = self.expand_one_record_from_atcc if self.atcc_database else self.expand_one_record_from_ncbi

    def run(self, gene_list_path: str):
        output_path = StrConverter.generate_result_file_name(gene_list_path, self.output_directory, 'gentamycin')
        gene_df = pd.read_csv(gene_list_path, sep='\t')
        gene_df[self.expand_headers] = gene_df.apply(lambda record: self.expand_one_record_func(record),
                                                     axis=1,
                                                     result_type='expand')
        gene_df.to_csv(output_path, sep='\t', index=None, header=True)
        prepare_consistency_file = StrConverter.generate_result_file_name(gene_list_path, self.output_directory,
                                                                          'gentamycin_consistency')
        consistency_df = self.generate_consistency_df(gene_df)
        consistency_df.to_csv(prepare_consistency_file, sep='\t', index=None, header=True)

    def generate_consistency_df(self, gene_df):
        if self.atcc_database:
            consistency_df = []
            for _, record in gene_df.iterrows():
                name = record[self.__header_gene__]
                sequences = record[self.__header_sequence__]
                if not sequences:
                    continue
                sequence_dict = json.loads(sequences)
                for tag, sequence in sequence_dict.items():
                    consistency_df.append({
                        'name': name + '-' + tag,
                        'gene': sequence
                    })
            consistency_df = pd.DataFrame(consistency_df)
        else:
            consistency_df = gene_df[['Gene', 'sequence']].rename({'Gene': 'name', 'sequence': 'sequence'})
        return consistency_df

    def expand_one_record_from_atcc(self, record: pd.Series):
        """
        :param record:
        :return: left gene, right gene, hit gene, sequence
        """
        if not record[self.__header_gene__].startswith('DR'):
            return tuple(['' for _ in self.expand_headers])
        else:

            left, right, direction = self.get_position(record['Locus'].strip())
            result = self.gene_position_helper.get_nearby_gene_based_by_range(left, right, direction)
            result = tuple([result.get(header) for header in self.expand_headers])
            return result

    def expand_one_record_from_ncbi(self, record: pd.Series):
        """
        :param record:
        :return: sequence
        """
        left, right, direction = self.get_position(record['Locus'].strip())
        sequence = self.ncbi_database.dna_code[left - 1:right]
        if direction == '-':
            sequence = get_opposite_dna(sequence[::-1])
        return sequence,

    @staticmethod
    def get_position(locus):
        matched = re.findall(r'(.+):(\d+)-(\d+)\((.)\)', locus)
        assert matched and len(matched) == 1, locus
        matched = matched[0]
        left, right, direction = int(matched[1]), int(matched[2]), matched[3]
        return left, right, direction
