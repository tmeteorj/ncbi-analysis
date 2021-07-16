from utils.gene_database import GeneDatabase
from utils.gene_util import get_opposite_dna


class GenePositionHelper:
    def __init__(self, database: GeneDatabase):
        self.database = database

    def get_nearby_gene_based_by_range(self, left, right, direction):
        result = {'related': [], 'hit': None, 'sequence': {'hit': '', 'related': []}}
        left_ge_id = self.database.find_first_greater_equal(left)
        right_ge_id = self.database.find_first_greater_equal(right)
        right_lt_id = right_ge_id - 1
        if left_ge_id == right_lt_id:
            result['hit'] = self.database.gene_segments[left_ge_id].gene
            result['sequence']['hit'] = self.database.get_sequence(left_ge_id)
        elif left_ge_id < right_lt_id:
            for idx in range(left_ge_id, right_lt_id + 1):
                result['related'].append(self.database.gene_segments[idx].gene)
                result['sequence']['related'].append(self.database.get_sequence(idx))
        else:
            if left_ge_id - 1 != right_lt_id:
                raise ValueError('left_ge_id-1!=right_lt_id')
            for idx in [right_lt_id, left_ge_id]:
                result['related'].append(self.database.gene_segments[idx].gene)
                result['sequence']['related'].append(self.database.get_sequence(idx))
        if direction == '-':
            self.change_to_opposite(result['sequence'])
        self.drop_empty_result(result)
        return result

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
