from enum import Enum

from utils.gene_database import GeneDatabase, GeneSegment
from utils.gene_util import get_opposite_dna


class RangeGenePositionType:
    Hit = 'hit'
    Homology = 'homology'
    Include = 'include'
    Related = 'related'


class GenePositionHelper:
    def __init__(self, database: GeneDatabase):
        self.database = database

    def get_nearby_gene_based_by_range(self, left, right, direction):
        left_ge_id = self.database.find_first_greater_equal(left)
        right_ge_id = self.database.find_first_greater_equal(right)
        right_lt_id = right_ge_id - 1
        if left_ge_id == right_lt_id:
            yield self.generate_result(left, right, left_ge_id, direction)
        elif left_ge_id < right_lt_id:
            covered_gene_idx = self.find_gene_cover_range(range(left_ge_id, right_lt_id + 1),
                                                          left,
                                                          right)
            if covered_gene_idx:
                yield self.generate_result(left, right, covered_gene_idx, direction)
            else:
                for idx in range(left_ge_id, right_lt_id + 1):
                    yield self.generate_result(left, right, idx, direction)
        else:
            if left_ge_id - 1 != right_lt_id:
                raise ValueError('left_ge_id-1!=right_lt_id')
            covered_gene_idx = self.find_gene_cover_range([right_lt_id, left_ge_id],
                                                          left,
                                                          right)
            if covered_gene_idx:
                yield self.generate_result(left, right, covered_gene_idx, direction)
            else:
                for idx in [right_lt_id, left_ge_id]:
                    yield self.generate_result(left, right, idx, direction)

    def generate_result(self, range_left, range_right, idx, direction):
        segment = self.database.gene_segments[idx]
        sequence = self.database.get_sequence(idx)
        if direction == '-':
            sequence = get_opposite_dna(sequence)
        return {
            'type': self.get_gene_range_type(range_left,
                                             range_right,
                                             segment.left,
                                             segment.right),
            'gene': segment.gene,
            'gene_left': segment.left,
            'gene_right': segment.right,
            'sequence': sequence
        }

    def find_gene_cover_range(self, idx_range, left, right):
        for idx in idx_range:
            gene_segment = self.database.gene_segments[idx]
            if gene_segment.left <= left and gene_segment.right >= right:
                return idx
        return None

    def drop_empty_result(self, result):
        if isinstance(result, dict):
            items = list(result.items())
            for k, v in items:
                if not v:
                    result.pop(k)
                elif isinstance(v, dict):
                    self.drop_empty_result(v)

    def change_to_opposite(self, result):
        if isinstance(result, dict):
            for tag, sequence in result.items():
                if isinstance(sequence, str):
                    result[tag] = get_opposite_dna(sequence[::-1])
                else:
                    self.change_to_opposite(result[tag])
        elif isinstance(result, list):
            for idx, sequence in enumerate(result):
                if isinstance(sequence, str):
                    result[idx] = get_opposite_dna(sequence[::-1])
                else:
                    self.change_to_opposite(result[idx])
        else:
            raise ValueError('change_to_opposite(%s)' % (type(result)))

    def get_gene_range_type(self, range_left, range_right, gene_left, gene_right):
        if range_left == gene_left and range_right == gene_right:
            return RangeGenePositionType.Hit
        overlap_size = self.get_overlap_size(range_left, range_right, gene_left, gene_right)
        overlap_rate = overlap_size * 100.0 / (range_right - range_left + 1)
        if overlap_rate >= 90.0:
            return RangeGenePositionType.Homology
        if gene_left <= range_left <= range_right <= gene_right:
            return RangeGenePositionType.Include
        else:
            return RangeGenePositionType.Related

    def get_diff_size(self, range_left, range_right, gene_left, gene_right):
        range_size = range_right - range_left + 1
        gene_size = gene_right - gene_left + 1
        diff_size = range_size - gene_size
        return diff_size if diff_size > 0 else -diff_size

    def get_overlap_size(self, range_left, range_right, gene_left, gene_right):
        if range_right < gene_left or range_left > gene_right:
            return 0
        if range_right < gene_right:
            return range_right - gene_left + 1
        else:
            return gene_right - range_left + 1
