import os
from dataclasses import dataclass

from utils.factories.logger_factory import LoggerFactory
from utils.ncbi_database import NCBIDatabase
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter


@dataclass
class GeneStreamAnalysis:
    data_path: str
    input_path: str
    output_directory: str
    mode: str = 'rna'
    limit: int = 200

    def __post_init__(self):
        self.inter_path = self.input_path if self.mode == 'inter' else None
        self.rna_path = self.input_path if self.mode == 'rna' else None
        file_name = os.path.basename(self.input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        suffix = 'stream_%d' % self.limit if self.mode == 'rna' else 'gene'
        self.result_path = os.path.join(self.output_directory, '%s_%s_result.txt' % (file_prefix, suffix))
        self.gene_reader = NCBIDatabase(self.data_path)
        self.logger = LoggerFactory()
        self.headers = {}
        self.inv_headers = []

    def get_utr_between(self, first, second):
        left = self.gene_reader.gene_segments[first].cds[1]
        right = self.gene_reader.gene_segments[second].cds[0] - 1
        return self.gene_reader.dna_code[left:right]

    def work_for_gene_index(self, index, start, end):
        gene_segment = self.gene_reader.gene_segments[index]
        assert gene_segment.cds[0] == min(start, end)
        assert gene_segment.cds[1] == max(start, end)
        seq = self.gene_reader.dna_code[gene_segment.cds[0] - 1:gene_segment.cds[1]]
        upstream = self.gene_reader.dna_code[max(gene_segment.cds[0] - self.limit - 1, 0):gene_segment.cds[0] - 1]
        downstream = self.gene_reader.dna_code[gene_segment.cds[1]:gene_segment.cds[1] + self.limit]
        if start > end:
            seq = get_opposite_dna(seq[::-1])
            upstream, downstream = get_opposite_dna(downstream[::-1]), get_opposite_dna(upstream[::-1])
        return seq, upstream, downstream

    def work_for_gene(self, gene_idx, gene_name, start, end, fw):
        if gene_name.find('->') >= 0:
            gene_name = gene_name[:gene_name.index('->')]
        if gene_name not in self.gene_reader.gene_name_segment_map:
            self.logger.info("%s not found in data" % gene_name)
            return
        cnt = 1
        fw.write('%d. %s\n' % (gene_idx, gene_name))
        for idx in self.gene_reader.gene_name_segment_map[gene_name]:
            seq, up, down = self.work_for_gene_index(idx, start, end)
            fw.write('%d)\n' % cnt)
            fw.write('position\t%d %s %d\n' % (
                self.gene_reader.gene_segments[idx].cds[0], '->' if start < end else '<-',
                self.gene_reader.gene_segments[idx].cds[1]))
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

    def generate_header(self, items):
        for idx, col_name in enumerate(items.strip().split('\t')):
            self.headers[col_name] = idx
            self.inv_headers.append(col_name)

    def run(self):
        with open(self.result_path, 'w', encoding='utf8') as fw:
            if self.mode == 'rna':
                lines = open(self.rna_path, 'r', encoding='utf8').readlines()
                self.generate_header(lines[0])
                for gene_idx, line in enumerate(lines[1:]):
                    items = line.split('\t')
                    gene_name, start, end = items[self.headers['gene']], int(
                        items[self.headers['map_start_pos']]), int(items[self.headers['map_end_pos']])
                    self.work_for_gene(gene_idx, gene_name.strip(), start, end, fw)
            elif self.mode == 'inter':
                self.check_inter(fw)
            else:
                raise ValueError(self.mode)
