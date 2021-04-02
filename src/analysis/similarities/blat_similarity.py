from analysis.similarities.base_similarity import BaseSimilarity


class BlatSimilarity(BaseSimilarity):
    def __init__(self, mid_limit: int = 10, end_limit: int = 2):
        self.end_limit = end_limit
        self.mid_limit = mid_limit

    def get_similairty(self, gene: str, database: str, offset: int):
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
