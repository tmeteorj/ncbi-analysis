import os

from utils.gene_file_util import GeneFileReader
from utils.str_util import StrConverter


class GeneRangeExtract:
    def __init__(self, data_path, range_path, output_directory):
        self.data_path = data_path
        self.range_path = range_path
        self.output_directory = output_directory
        file_name = os.path.basename(range_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory, '%s_range_result.txt' % file_prefix)
        self.gene_reader = GeneFileReader(self.data_path)

    def generate_header(self, items):
        for idx, col_name in enumerate(items.strip().split('\t')):
            self.headers[col_name] = idx
            self.inv_headers.append(col_name)

    def run(self):
        self.gene_reader.build_information()
        dna_code = self.gene_reader.dna_code
        with open(self.result_path, 'w', encoding='utf8') as fw:
            self.extract_range(dna_code, fw)

    def extract_range(self, dna_code, fw):
        lines = [line.strip() for line in open(self.range_path, 'r', encoding='utf8')]
        last_end = 0
        gene_idx = 0
        fw.write('name\trange\n')
        for line in lines[1:]:
            items = line.split('\t')
            locus = items[0]
            locus = locus.split(':')[1]
            left, right = map(int, locus.split('-'))
            if last_end < left - 1:
                gene_idx += 1
                fw.write('gene_%d\tNC_000913.3:%d-%d\n' % (gene_idx, last_end, left - 1))
            last_end = right
