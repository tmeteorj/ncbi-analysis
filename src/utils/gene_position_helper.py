from utils.gene_database import GeneDatabase


class GenePositionHelper:
    def __init__(self, database: GeneDatabase):
        self.database = database

    def get_nearby_gene_based_by_range(self, left, right):
        result = {'left': None, 'right': None, 'hit': None, 'sequence': {}}
        left_ge_id = self.database.find_first_greater_equal(left)
        right_ge_id = self.database.find_first_greater_equal(right)
        right_lt_id = right_ge_id - 1
        if left_ge_id == right_lt_id:
            result['hit'] = self.database.gene_segments[left_ge_id].gene
            result['sequence']['hit'] = self.database.get_sequence(left_ge_id)
        elif left_ge_id < right_lt_id:
            assert left_ge_id + 1 == right_lt_id
            result['left'] = self.database.gene_segments[left_ge_id].gene
            result['right'] = self.database.gene_segments[right_lt_id].gene
            result['sequence']['left'] = self.database.get_sequence(left_ge_id)
            result['sequence']['right'] = self.database.get_sequence(right_lt_id)
        else:
            assert left_ge_id - 1 == right_lt_id
            result['left'] = self.database.gene_segments[right_lt_id].gene
            result['right'] = self.database.gene_segments[left_ge_id].gene
            result['sequence']['left'] = self.database.get_sequence(right_lt_id)
            result['sequence']['right'] = self.database.get_sequence(left_ge_id)
        return result
