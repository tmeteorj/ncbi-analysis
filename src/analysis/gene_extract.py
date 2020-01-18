import os
import traceback

from utils.gene_file_util import GeneFileReader
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter


class GeneExtract:
    def __init__(self, data_path, rna_path, output_directory, gene_extract_based='gene', left_idx=-2, right_idx=-1):
        self.data_path = data_path
        self.rna_path = rna_path
        self.output_directory = output_directory
        file_name = os.path.basename(rna_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory, '%s_extract_result.txt' % file_prefix)
        self.gene_extract_based = gene_extract_based
        self.gene_reader = GeneFileReader(self.data_path)
        self.left_idx = left_idx
        self.right_idx = right_idx
        self.headers = {}
        self.inv_headers = []

    def generate_header(self, items):
        for idx, col_name in enumerate(items.strip().split('\t')):
            self.headers[col_name] = idx
            self.inv_headers.append(col_name)

    def run(self):
        self.gene_reader.build_information()
        dna_code = self.gene_reader.dna_code
        with open(self.result_path, 'w', encoding='utf8') as fw:
            if self.gene_extract_based == 'gene':
                self.extract_sequence_based_on_gene(dna_code, fw)
            elif self.gene_extract_based == 'range':
                self.extract_sequence_based_on_range(dna_code, fw)

    def extract_sequence_based_on_gene(self, dna_code, fw):
        fw.write('No\tgene\tfrom\t\tend\tproduct\tsequence\n')
        for gene_idx, gene in enumerate(open(self.rna_path)):
            gene = gene.strip()
            succ = False
            for idx in self.gene_reader.gene_name_segment_map.get(gene, []):
                gene_segment = self.gene_reader.gene_segments[idx]
                succ = True
                start = gene_segment.cds[0]
                end = gene_segment.cds[1]
                product = gene_segment.product
                sequence = dna_code[start - 1:end]
                fw.write('d%d\t%s\t%s\t%s\t%s\t%s\n' % (
                    gene_idx + 1, gene, start, end, product, sequence))
            if not succ:
                print('%s not found in %s' % (gene, self.data_path))

    def extract_sequence_based_on_range(self, dna_code, fw):
        lines = [line.strip() for line in open(self.rna_path, 'r', encoding='utf8')]
        self.generate_header(lines[0])
        fw.write(lines[0] + '\n')
        for line in lines[1:]:
            result = {}
            infos = line.strip().split('\t')
            for idx, info in enumerate(infos):
                result[self.inv_headers[idx]] = info
            if result.get('sequence', '') == '':
                try:
                    a, b = map(int, [infos[self.left_idx], infos[self.right_idx]])
                    left = min(a, b)
                    right = max(a, b)
                    direction = a < b
                    # id start from 0
                    left -= 1
                    right -= 1
                    if not direction:
                        left += 1
                        right += 1
                    dna = dna_code[left:right]
                    if not direction:
                        result['sequence'] = get_opposite_dna(dna[::-1])
                    else:
                        result['sequence'] = dna
                except:
                    print(infos)
                    traceback.print_exc()
            fw.write(self.extract_output(result) + '\n')

    def extract_output(self, result):
        output = []
        for name in self.inv_headers:
            output.append(result.get(name, ''))
        return '\t'.join(output)
