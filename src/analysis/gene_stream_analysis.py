import os

from utils.factories.logger_factory import LoggerFactory
from utils.gene_file_util import GeneFileReader
from utils.str_util import StrConverter


class GeneStreamAnalysis:
    def __init__(self, data_path, input_path, output_directory, mode='rna'):
        self.mode = mode
        self.data_path = data_path
        self.inter_path = input_path if mode == 'inter' else None
        self.rna_path = input_path if mode == 'rna' else None
        self.output_directory = output_directory
        file_name = os.path.basename(input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        suffix = 'stream' if mode == 'rna' else 'gene'
        self.result_path = os.path.join(self.output_directory, '%s_%s_result.txt' % (file_prefix, suffix))
        self.gene_reader = GeneFileReader(self.data_path)
        self.logger = LoggerFactory()

    def get_utr_between(self, first, second):
        left = self.gene_reader.gene_segments[first].cds[1]
        right = self.gene_reader.gene_segments[second].cds[0] - 1
        return self.gene_reader.dna_code[left:right]

    def work_for_gene_index(self, index):
        gene_segment = self.gene_reader.gene_segments[index]
        seq = self.gene_reader.dna_code[gene_segment.cds[0] - 1:gene_segment.cds[1]]
        upstream = self.get_utr_between(index - 1, index) if index > 0 else None
        downstream = self.get_utr_between(index, index + 1) if index < len(self.gene_reader.gene_segments) - 1 else None
        return seq, upstream, downstream

    def work_for_gene(self, gene_idx, gene_name, fw):
        if gene_name not in self.gene_reader.gene_name_segment_map:
            self.logger.info("%s not found in data" % gene_name)
            return
        cnt = 1
        fw.write('%d. %s\n' % (gene_idx, gene_name))
        for idx in self.gene_reader.gene_name_segment_map[gene_name]:
            seq, up, down = self.work_for_gene_index(idx)
            fw.write('%d)\n' % cnt)
            fw.write('position\t%s\n' % ('-'.join(map(str, self.gene_reader.gene_segments[idx].cds))))
            fw.write('product\t%s\n' % self.gene_reader.gene_segments[idx].product)
            fw.write('GeneID\t%s\n' % self.gene_reader.gene_segments[idx].gene_id)
            fw.write('stream\t%s\n' % seq)
            if up: fw.write('upstream\t%s\n' % up)
            if down: fw.write('downstream\t%s\n' % down)
            fw.write('\n')
            cnt += 1

    def check_inter(self, fw):
        for line in open(self.inter_path, 'r', encoding='utf8'):
            line = line.strip()
            if line == '': continue
            left, right = map(int, line.split(','))
            up, down = None, None
            for gene_segment in self.gene_reader.gene_segments:
                if max(gene_segment.cds) < left:
                    if not up or max(up.cds) < max(gene_segment.cds):
                        up = gene_segment
                if min(gene_segment.cds) > right:
                    if not down or min(down.cds) > min(gene_segment.cds):
                        down = gene_segment
            fw.write('%s:\n' % line)
            if up:
                fw.write('up-gene\t%s\nup-position\t%s\nup-product\t%s\n' % (
                    up.gene, '-'.join(map(str, up.cds)), up.product))
            if down:
                fw.write('down-gene\t%s\ndown-position\t%s\ndown-product\t%s\n' % (
                    down.gene, '-'.join(map(str, down.cds)), down.product))
            fw.write('\n')

    def run(self):
        self.gene_reader.build_information()
        with open(self.result_path, 'w', encoding='utf8') as fw:
            if self.mode == 'rna':
                for gene_idx, gene_name in enumerate(open(self.rna_path, 'r', encoding='utf8')):
                    self.work_for_gene(gene_idx, gene_name.strip(), fw)
            elif self.mode == 'inter':
                self.check_inter(fw)
            else:
                raise ValueError(self.mode)
