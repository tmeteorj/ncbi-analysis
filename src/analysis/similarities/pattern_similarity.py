import re

from analysis.similarities.base_similarity import BaseSimilarity


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


class PatternSimilarity(BaseSimilarity):
    def __init__(self, match_pattern: MatchPattern):
        self.match_pattern = match_pattern

    def get_similarity(self, gene: str, database: str, offset: int):
        if self.match_pattern is None:
            return 0
        gene_len = len(gene)
        target_gene = database[offset:offset + gene_len]
        if not re.match(self.match_pattern.must_pattern, target_gene):
            return 0
        score = self.match_pattern.must_score
        for optional_pattern in self.match_pattern.option_patterns:
            if re.match(optional_pattern[0], target_gene):
                score += optional_pattern[1]
        return score
