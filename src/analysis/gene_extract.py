import os

from utils.gene_file_util import GeneFileReader
from utils.str_util import StrConverter


class GeneExtract:
    def __init__(self, data_path, rna_path, output_directory, gene_extract_based='gene'):
        self.data_path = data_path
        self.rna_path = rna_path
        self.output_directory = output_directory
        file_name = os.path.basename(rna_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory, '%s_extract_result.txt' % file_prefix)
        self.gene_extract_based = gene_extract_based
        self.gene_reader = GeneFileReader(self.data_path)

    def run(self):
        self.gene_reader.build_information()
        dna_code = self.gene_reader.dna_code
        with open(self.result_path, 'w') as fw:
            if self.gene_extract_based == 'gene':
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
            elif self.gene_extract_based == 'range':
                for line in open(self.rna_path):
                    fw.write(line.strip('\n'))
                    info = line.strip().split('\t')
                    if len(info) == 5:
                        a, b = map(int, info[-2:])
                        left = min(a, b)
                        right = max(a, b)
                        direction = a < b
                        if direction:
                            right += 2
                        else:
                            left -= 2
                        dna = dna_code[left:right + 1]
                        if not direction:
                            fw.write('\t' + dna[::-1])
                        else:
                            fw.write('\t' + dna)
                    fw.write('\n')
