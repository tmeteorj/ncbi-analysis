import re
from typing import Tuple

from analysis.models.match_pattern import MatchPattern
from analysis.similarities.base_similarity import BaseSimilarity


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

    def rendering_sequence(self, gene: str, database: str, offset: int) -> Tuple[list, list, list]:
        sequence_gene = []
        sequence_target = []
        sequence = []
        tot = len(gene)
        for i in range(tot):
            sequence_gene.append(gene[i])
            sequence_target.append(database[i + offset])
            if not self.should_change(gene[i], database[i + offset]):
                sequence.append('*')
            else:
                sequence.append('.')
        return sequence_gene, sequence_target, sequence
