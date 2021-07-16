import re
import traceback
from enum import Enum

from experiment_config import *
from utils.factories.logger_factory import LoggerFactory
from utils.gene_database import GeneSegment, GeneDatabase


class GeneDataPartType(Enum):
    HeaderPart = 0
    GeneSegmentPart = 1
    DNAPart = 2


class GeneDataLineType(Enum):
    SourceLine = 0
    GeneSegmentStart = 1
    DNAStart = 2
    DNAEnd = 3
    Other = 4


class NCBIGeneSegment(GeneSegment):
    def __init__(self):
        self.xref = {}
        self.cds = None
        self.gene_id = None

        self.product = None
        self.gene = None
        self.protein_id = None
        self.codon_start = None
        self.transl_table = None
        self.gene_synonym = None
        self.locus_tag = None
        self.translation = None

    def extract_attribute(self, line):
        for attr in ['product', 'gene', 'protein_id', 'codon_start', 'transl_table', 'gene_synonym', 'locus_tag',
                     'translation']:
            if line.startswith('/' + attr + '='):
                self.__dict__[attr] = line[len(attr) + 2:].strip('"')
        if line.startswith("/db_xref="):
            key, value = line.lstrip("/db_xref=").strip("\"").split(':')
            if key.lower() == 'geneid':
                value = value.split("\"")[0]
                self.gene_id = int(re.sub(r'[^0-9]', '', value.lower()))
            else:
                self.xref[key] = value

    def __str__(self):
        return '%s-%s\t%s' % (
            self.cds[0],
            self.cds[1],
            ExperimentConfig.VALUE_UNKNOWN if self.product is None else self.product)


class NCBIDatabase(GeneDatabase):
    def __init__(self, file_path, ignore_gene=False, enable_debug_info=False):
        self.ignore_gene = ignore_gene
        self.gene_segments = []
        self.dna_code = []
        self.gene_name_segment_map = {}
        self.source = None

        self.enable_debug_info = enable_debug_info
        self.file_path = file_path
        self.logger = LoggerFactory(1)

    def initialize(self):
        part_status = GeneDataPartType.HeaderPart
        data = []
        line_type = None
        for line_index, line in enumerate(open(self.file_path, 'r', encoding='utf8')):
            line_type = self.check_line_type(line, part_status)
            if line_type == GeneDataLineType.SourceLine:
                self.source = ' '.join(re.split(r'\s+', line)[1:])
            elif line_type == GeneDataLineType.GeneSegmentStart:
                part_status = GeneDataPartType.GeneSegmentPart
                self.parse_gene_segment(data)
            elif line_type == GeneDataLineType.DNAStart:
                part_status = GeneDataPartType.DNAPart
                self.parse_gene_segment(data)
            elif line_type == GeneDataLineType.DNAEnd:
                break

            if part_status == GeneDataPartType.GeneSegmentPart:
                data.append(line)
            elif part_status == GeneDataPartType.DNAPart and line_type == GeneDataLineType.Other:
                items = re.split(r'\s+', line.strip())
                for val in items[1:]:
                    self.dna_code.append(val)
                self.logger.time_distance = 10
            if self.enable_debug_info:
                self.logger.info_per_time("LineNo = %d, Added Gene Num = %d, Last Sample = %s" % (
                    line_index, len(self.gene_segments),
                    self.gene_segments[-1].__dict__ if len(self.gene_segments) > 0 else ""))
        if part_status != GeneDataPartType.DNAPart and line_type != GeneDataLineType.DNAEnd:
            return False
        self.dna_code = ''.join(self.dna_code)
        check_order = None
        warning_num = 0
        for idx, gene_segment in enumerate(self.gene_segments):
            name = gene_segment.gene
            if check_order is not None and check_order > min(gene_segment.cds):
                warning_num += 1
            check_order = max(gene_segment.cds)
            if name not in self.gene_name_segment_map:
                self.gene_name_segment_map[name] = []
            self.gene_name_segment_map[name].append(idx)
        self.logger.info("Total Gene Segment Number = %d, Total Gene Name Count = %d" % (
            len(self.gene_segments), len(self.gene_name_segment_map)))
        return True

    def parse_gene_segment(self, data):
        if data is None or len(data) == 0 or self.ignore_gene:
            return
        gene_segment = NCBIGeneSegment()
        last_line = ''
        success = True
        complement = None
        for line in data:
            try:
                line_type = self.check_line_type(line, GeneDataPartType.GeneSegmentPart)
                line = line.strip()
                if line_type == GeneDataLineType.GeneSegmentStart:
                    tag, complement = re.split(r'\s+', line)
                    inter = list(map(lambda arg: int(arg.strip('<>')),
                                     complement.lstrip('complement(').rstrip(')').split('..')))
                    gene_segment.cds = inter
                    assert (inter[0] < inter[1])
                else:
                    if line[0] == '/':
                        last_line = line
                    else:
                        last_line += ' ' + line
                    gene_segment.extract_attribute(last_line)
            except:
                self.logger.info(line)
                if not complement or (
                        not complement.startswith('join') and not complement.startswith('complement(join')):
                    traceback.print_exc()
                success = False
                break
        if success:
            self.gene_segments.append(gene_segment)
            gene_segment.left, gene_segment.right = gene_segment.cds[0], gene_segment.cds[1]
        data.clear()

    @staticmethod
    def check_line_type(line: str, part_status):
        strip_line = line.strip()
        if part_status == GeneDataPartType.HeaderPart:
            if strip_line.startswith(ExperimentConfig.VALUE_SOURCE_START):
                return GeneDataLineType.SourceLine
            elif strip_line.startswith(ExperimentConfig.VALUE_GENE_START) or strip_line.startswith(
                    ExperimentConfig.VALUE_REPEAT_REGION_START):
                return GeneDataLineType.GeneSegmentStart
        elif part_status == GeneDataPartType.GeneSegmentPart:
            if strip_line.startswith(ExperimentConfig.VALUE_GENE_START) or strip_line.startswith(
                    ExperimentConfig.VALUE_REPEAT_REGION_START):
                return GeneDataLineType.GeneSegmentStart
            elif line[0] != ' ':
                return GeneDataLineType.DNAStart
        elif part_status == GeneDataPartType.DNAPart:
            if strip_line.startswith(ExperimentConfig.VALUE_DNA_PART_END):
                return GeneDataLineType.DNAEnd
        return GeneDataLineType.Other


if __name__ == '__main__':

    file_path = 'D:/Workspace/ncbi-analysis/data/rna_analysis/rna_download_data/NC_000913.3.txt'
    gene = NCBIDatabase(file_path)
    gene.initialize()
    b_dict = {}
    for gene_segment in gene.gene_segments:
        if gene_segment.gene_id:
            b_dict[gene_segment.locus_tag] = [gene_segment.gene, gene_segment.gene_id]
    heads = None
    with open('D:/Workspace/ncbi-analysis/R/Path_Information_Sample.txt', 'w', encoding='utf8') as fw:
        fw.write('PathId\tPath Name\tLocusTag\tGene\tGeneId\n')
        for line in open('D:/Workspace/ncbi-analysis/R/ID_Eco_Path.txt', 'r', encoding='utf8'):
            items = line.strip().split('\t')
            if heads is None:
                heads = items
                # keggid locus path
            elif items[1] in ['b2397',
                              'b2815',
                              'b2590',
                              'b0745',
                              'b1231',
                              'b0216',
                              'b2402',
                              'b3171',
                              'b1975',
                              'b4370',
                              'b0971',
                              'b3798',
                              'b3761',
                              'b1911',
                              'b1977',
                              'b0749']:
                items[1], items[2] = items[2], items[1]
                items.extend(b_dict.get(items[2], ['', '']))
                items = map(str, items)
                fw.write('\t'.join(items) + '\n')
