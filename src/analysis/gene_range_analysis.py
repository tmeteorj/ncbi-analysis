import os

from utils.gene_file_util import GeneFileReader
from utils.str_util import StrConverter


class GeneRangeExtract:
    def __init__(self, data_path, output_directory):
        self.data_path = data_path
        self.output_directory = output_directory
        file_name = os.path.basename(data_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory, '%s_range_result.txt' % file_prefix)
        self.gene_reader = GeneFileReader(self.data_path)

    def generate_header(self, items):
        for idx, col_name in enumerate(items.strip().split('\t')):
            self.headers[col_name] = idx
            self.inv_headers.append(col_name)

    def run(self):
        self.gene_reader.build_information()
        with open(self.result_path, 'w', encoding='utf8') as fw:
            last_end = 0
            gene_idx = 0
            fw.write('name\trange\n')
            for gene_segment in self.gene_reader.gene_segments:
                left, right = gene_segment.cds
                if last_end < left - 1:
                    gene_idx += 1
                    fw.write('region_%d\tNC_000913.3:%d-%d\n' % (gene_idx, last_end + 1, left - 1))
                if left > right:
                    print(gene_segment)
                last_end = right
