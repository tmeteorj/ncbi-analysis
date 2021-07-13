import re

import pandas as pd

from utils.factories.logger_factory import LoggerFactory

logger = LoggerFactory()


class ColiGeneSegment:
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


class ColiDatabase(object):
    def __init__(self, coli_path):
        self.segments = []
        buff = []
        for line in open(coli_path, 'r'):
            if line.startswith('>lcl'):
                if len(buff) > 0:
                    self.segments.append(ColiGeneSegment(buff))
                    buff.clear()
            buff.append(line.strip())
        if len(buff) > 0:
            self.segments.append(ColiGeneSegment(buff))
        self.segments.sort(key=lambda arg: arg.start)
        logger.info('Segments Count = %d' % len(self.segments))

    def find_first_greater_equal(self, pos):
        start, end = 0, len(self.segments) - 1
        while start < end:
            mid = (start + end) // 2
            if self.segments[mid].start < pos:
                start = mid + 1
            elif self.segments[mid].start >= pos:
                end = mid
        if self.segments[end].start >= pos:
            return end
        else:
            return end + 1
