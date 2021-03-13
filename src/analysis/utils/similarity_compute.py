import re


class MatchPattern:
    must_pattern: str = None
    option_patterns: list = None
    must_score: int = 0

    def __init__(self, rna, conditions):
        must = conditions['must']
        self.must_pattern, self.must_score = self.generate_pattern(rna, must)
        self.option_patterns = []
        for optional in conditions['optional']:
            optional = [optional]
            optional.extend(must)
            optional_pattern, optioanl_score = self.generate_pattern(rna, optional)
            optioanl_score -= self.must_score
            self.option_patterns.append((optional_pattern, optioanl_score))

    def generate_pattern(self, rna, conditions):
        rna_len = len(rna)
        conditions.sort(key=lambda arg: arg['offset'] if arg['offset'] >= 0 else rna_len + arg['offset'])
        gen_pattern = ''
        score = 0
        index = 0
        for condition in conditions:
            offset, length = condition['offset'], condition['length']
            if offset < 0:
                offset = rna_len + offset
            if offset == 0:
                gen_pattern += '^'
            if offset > index:
                gen_pattern += '.+'
            gen_pattern += self.update_regex(rna[offset:offset + length])
            index = offset + length
            if index == rna_len:
                gen_pattern += '$'
            score += length
        if index != rna_len:
            gen_pattern += '.+'
        return gen_pattern, score

    def update_regex(self, pattern: str):
        up_pattern = ''
        pattern = pattern.lower()
        for c in pattern:
            if c == 'c':
                up_pattern += '(c|t)'
            else:
                up_pattern += c
        return up_pattern


def compute_text_distance_similarity(gene: str, database: str, offset: int, continuous_mismatch_limit: int = None):
    tot = len(gene)
    dp = [[99999 for _ in range(tot + 1)] for _ in range(tot + 1)]
    dp[0][0] = 0
    for i in range(1, tot + 1):
        gene_a = gene[i - 1]
        min_step = 1000000
        for j in range(1, tot + 1):
            gene_b = database[offset + j - 1]
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + should_change(gene_a, gene_b))
            min_step = min(dp[i][j] + abs(i - j), min_step)
    score = float(tot - dp[tot][tot])
    if continuous_mismatch_limit is not None:
        i, j = tot, tot
        mismatch = 0
        while i > 0 or j > 0:
            gene_a, gene_b = gene[i - 1] if i > 0 else '.', database[j + offset - 1] if j > 0 else '.'
            if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + should_change(gene[i - 1],
                                                                                database[j + offset - 1]):

                if should_change(gene[i - 1], database[j + offset - 1]) != 0:
                    mismatch += 1
                else:
                    mismatch = 0
                i, j = i - 1, j - 1
            elif dp[i][j] == dp[i - 1][j] + 1:
                mismatch += 1
                i -= 1
            elif dp[i][j] == dp[i][j - 1] + 1:
                mismatch += 1
                j -= 1
            else:
                raise ValueError('Should not go here!')
            if mismatch >= continuous_mismatch_limit:
                return 0, dp
    return score, dp


def compute_consistency_similarity(gene: str, database: str, offset: int, max_patience: int):
    tot = len(gene)
    score = 0.0
    same = 0
    cur_score = 0
    score_queue = []
    for i in range(tot):
        if not should_change(gene[i], database[i + offset]):
            cur_score += 1
            same += 1
            if i == tot - 1:
                score_queue.append([cur_score, tot])
        else:
            score_queue.append([cur_score, i])
            cur_score = 0
        score = max(score, cur_score)
    score_merge_idx = [-1, -1]
    for idx in range(len(score_queue)):
        left = score_queue[idx][1] - score_queue[idx][0]
        total_score = 0
        for width in range(max_patience + 1):
            if width + idx < len(score_queue):
                total_len = score_queue[idx + width][1] - left
                total_score += score_queue[idx + width][0]
                if total_len - total_score > max_patience:
                    break
                if score < total_score:
                    score = total_score
                    score_merge_idx = [idx, idx + width]
    return score, score_queue, score_merge_idx


def compute_pattern_similarity(gene: str, database: str, offset: int, match_pattern: MatchPattern):
    if match_pattern is None:
        return 0
    gene_len = len(gene)
    target_gene = database[offset:offset + gene_len]
    if not re.match(match_pattern.must_pattern, target_gene):
        return 0
    score = match_pattern.must_score
    for optional_pattern in match_pattern.option_patterns:
        if re.match(optional_pattern[0], target_gene):
            score += optional_pattern[1]
    return score


def compute_blat_similarity(gene: str, database: str, offset: int):
    def search_dfs(pos_gene, pos_data, insert_data):
        if pos_gene < 4:
            sequence_matched_length = 1
            cond = False
            while pos_gene < 4 and pos_data < len(database):
                while should_change(gene[pos_gene], database[pos_data]) > 0:
                    sequence_matched_length = 0
                    insert_data += 1
                    pos_data += 1
                    if insert_data > 4 or pos_data >= len(database):
                        return False, None
                if sequence_matched_length > 0:
                    cond = True
                sequence_matched_length += 1
                pos_gene += 1
                pos_data += 1
            if not cond:
                return False, None
            flag, pos_data_end = search_dfs(4, pos_data + 1, 1)
            return flag, pos_data_end
        elif pos_gene == 4:
            if insert_data > 20 or pos_data >= len(database):
                return False, None
            while should_change(gene[pos_gene], database[pos_data]) > 0:
                pos_data += 1
                insert_data += 1
                if pos_data > len(database) or insert_data > 20:
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
                while should_change(gene[pos_gene], database[pos_data]) > 0:
                    sequence_matched_length = 0
                    insert_data += 1
                    pos_data += 1
                    if insert_data > 4 or pos_data >= len(database):
                        return False, None
                if sequence_matched_length > 0:
                    cond = True
                sequence_matched_length += 1
                pos_gene += 1
                pos_data += 1
            return cond, pos_data

    if should_change(gene[0], database[offset]) > 0:
        return 0, None
    flag, pos_data_end = search_dfs(1, offset + 1, 0)
    return flag, pos_data_end


def should_change(a, b):
    if a == b:
        return 0
    elif a == 'c' and b == 't':
        return 0
    return 1
