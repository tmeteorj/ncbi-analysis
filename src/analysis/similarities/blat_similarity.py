from typing import Tuple, Any

from analysis.similarities.base_similarity import BaseSimilarity


class BlatSimilarity(BaseSimilarity):

    def __init__(self, mid_limit: int = 10, end_limit: int = 2):
        self.end_limit = end_limit
        self.mid_limit = mid_limit

    def rendering_sequence(self, gene: str, database: str, offset: int) -> Tuple[list, list, list]:
        sequence_gene = []
        sequence_target = []
        sequence = []
        flag, pos_data_end = self.get_similarity(gene, database, offset)
        pos_data = offset
        pos_gene = 0
        while pos_gene < 4:
            if self.should_change(gene[pos_gene], database[pos_data]) > 0:
                sequence_gene.append('-')
                sequence_target.append(database[pos_data])
                sequence.append('.')
                pos_data += 1
            else:
                sequence_gene.append(gene[pos_gene])
                sequence_target.append(database[pos_data])
                sequence.append('*')
                pos_gene += 1
                pos_data += 1
        rev_pos_gene = 7
        rev_pos_data = pos_data_end - 1
        rev_sequence_gene = []
        rev_sequence_target = []
        rev_sequence = []
        while rev_pos_gene > 3:
            if self.should_change(gene[rev_pos_gene], database[rev_pos_data]) > 0:
                rev_sequence_gene.append('-')
                rev_sequence_target.append(database[rev_pos_data])
                rev_sequence.append('.')
                rev_pos_data -= 1
            else:
                rev_sequence_gene.append(gene[rev_pos_gene])
                rev_sequence_target.append(database[rev_pos_data])
                rev_sequence.append('*')
                rev_pos_gene -= 1
                rev_pos_data -= 1
        while pos_data <= rev_pos_data:
            sequence_gene.append('-')
            sequence_target.append(database[pos_data])
            sequence.append('.')
            pos_data += 1
        sequence_gene.extend(rev_sequence_gene[::-1])
        sequence_target.extend(rev_sequence_target[::-1])
        sequence.extend(rev_sequence[::-1])
        return sequence_gene, sequence_target, sequence

    def get_similarity(self, gene: str, database: str, offset: int) -> Tuple[float, Any]:
        end_limit = 2
        mid_limit = 10

        def search_dfs(pos_gene, pos_data, insert_data):
            if pos_gene < 4:
                sequence_matched_length = 1
                cond = False
                while pos_gene < 4 and pos_data < len(database):
                    while self.should_change(gene[pos_gene], database[pos_data]) > 0:
                        sequence_matched_length = 0
                        insert_data += 1
                        pos_data += 1
                        if insert_data > end_limit or pos_data >= len(database):
                            return False, None
                    if sequence_matched_length > 0:
                        cond = True
                    sequence_matched_length += 1
                    pos_gene += 1
                    pos_data += 1
                if not cond:
                    return False, None
                flag, pos_data_end = search_dfs(4, pos_data + 3, 3)
                return flag, pos_data_end
            elif pos_gene == 4:
                if insert_data > mid_limit or pos_data >= len(database):
                    return False, None
                while self.should_change(gene[pos_gene], database[pos_data]) > 0:
                    pos_data += 1
                    insert_data += 1
                    if pos_data >= len(database) or insert_data > mid_limit:
                        return False, None
                flag, pos_data_end = search_dfs(5, pos_data + 1, 0)
                if flag:
                    return flag, pos_data_end
                else:
                    flag, pos_data_end = search_dfs(4, pos_data + 1, insert_data + 1)
                    return flag, pos_data_end
            else:
                sequence_matched_length = 1
                cond = False
                while pos_gene < 8 and pos_data < len(database):
                    while self.should_change(gene[pos_gene], database[pos_data]) > 0:
                        sequence_matched_length = 0
                        insert_data += 1
                        pos_data += 1
                        if insert_data > end_limit or pos_data >= len(database):
                            return False, None
                    if sequence_matched_length > 0:
                        cond = True
                    sequence_matched_length += 1
                    pos_gene += 1
                    pos_data += 1
                return cond, pos_data

        if self.should_change(gene[0], database[offset]) > 0:
            return 0, None
        flag, pos_data_end = search_dfs(1, offset + 1, 0)
        return flag, pos_data_end
