import os
from dataclasses import dataclass

from utils.gene_file_util import GeneFileReader
from utils.str_util import StrConverter


@dataclass
class GeneRangeExtract:
    data_path: str
    output_directory: str

    def __post_init__(self):
        file_name = os.path.basename(self.data_path)
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
            region_idx = 0
            fw.write('name\trange\tlocus_tag\n')
            for gene_idx, gene_segment in enumerate(self.gene_reader.gene_segments):
                left, right = gene_segment.cds
                if last_end < left - 1:
                    region_idx += 1
                    fw.write('region_%d\t%d-%d\n' % (region_idx, last_end + 1, left - 1))
                fw.write('gene_%d\t%d-%d\t%s\n' % (gene_idx + 1, left, right, gene_segment.locus_tag))
                last_end = right
            total_len = len(self.gene_reader.dna_code)
            if last_end < total_len:
                region_idx += 1
                fw.write('region_%d\t%d-%d\n' % (region_idx, last_end + 1, total_len))
