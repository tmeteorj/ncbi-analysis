import re

import pandas as pd

from utils.factories.logger_factory import LoggerFactory
from utils.gene_database import GeneSegment, GeneDatabase

logger = LoggerFactory()


class ATCCGeneSegment(GeneSegment):
    location: str = None
    locus_tag: str = None
    start: int = None
    end: int = None
    gbkey: str = None
    gene: str = None
    sequence: str = None

    def __init__(self, buff):
        for item in buff[0].split():
            for attr in ['locus_tag', 'location', 'gbkey', 'gene']:
                matched = re.findall(rf'^\[{attr}=(.+)\]$', item, re.IGNORECASE)
                if matched:
                    assert len(matched) == 1
                    setattr(self, attr, matched[0])
        self.parse_location()
        self.sequence = ''.join(buff[1:]).lower()
        self.gene = self.gene if self.gene else f'Unknown:{self.locus_tag}'

    def parse_location(self):
        if self.location:
            for pattern in [r'complement\((\d+)\.\.(\d+)\)', r'(\d+)\.\.(\d+)']:
                matched = re.findall(pattern, self.location, re.IGNORECASE)
                if matched:
                    assert len(matched) == 1
                    self.start = int(matched[0][0])
                    self.end = int(matched[0][1])
                    break

    def to_dict(self):
        return {
            'locus_tag': self.locus_tag,
            'start': self.start,
            'end': self.end,
            'gbkey': self.gbkey,
            'gene': self.gene,
            'location': self.location,
            'sequence': self.sequence
        }


class ATCCDatabase(GeneDatabase):

    def __init__(self, coli_path):
        self.gene_segments: ATCCGeneSegment = []
        buff = []
        for line in open(coli_path, 'r'):
            if line.startswith('>lcl'):
                if len(buff) > 0:
                    self.gene_segments.append(ATCCGeneSegment(buff))
                    buff.clear()
            buff.append(line.strip())
        if len(buff) > 0:
            self.gene_segments.append(ATCCGeneSegment(buff))
        self.gene_segments.sort(key=lambda arg: arg.start)
        logger.info('Segments Count = %d' % len(self.gene_segments))

    def get_sequence(self, segment_id=None, left=None, right=None):
        return self.gene_segments[segment_id].sequence
