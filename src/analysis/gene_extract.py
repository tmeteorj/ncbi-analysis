from utils.gene_file_util import GeneFileReader


class GeneExtract:
    def __init__(self, data_path, rna_path, output_path, based='gene'):
        self.data_path = data_path
        self.rna_path = rna_path
        self.output_path = output_path
        self.based = based
        self.gene_reader = GeneFileReader(self.data_path)

    def run(self):
        self.gene_reader.build_information()
        dna_code = self.gene_reader.dna_code
        with open(self.output_path, 'w') as fw:
            if self.based == 'gene':
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
            elif self.based == 'internal':
                raise NotImplementedError()
