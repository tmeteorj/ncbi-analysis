from typing import Tuple, Any

from analysis.similarities.base_similarity import BaseSimilarity


class DirectMatchSimilarity(BaseSimilarity):

    def get_similarity(self, gene: str, database: str, offset: int) -> Tuple[float, Any]:
        score = 0.0
        tot = len(gene)
        for i in range(tot):
            score += self.should_change(gene[i], database[i + offset])
        score = tot - score
        return score, None
