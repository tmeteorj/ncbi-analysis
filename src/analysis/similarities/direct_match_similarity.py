from typing import Tuple, Any

from analysis.similarities.base_similarity import BaseSimilarity


class DirectMatchSimilarity(BaseSimilarity):

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

    def get_similarity(self, gene: str, database: str, offset: int) -> Tuple[float, Any]:
        score = 0.0
        tot = len(gene)
        for i in range(tot):
            score += self.should_change(gene[i], database[i + offset])
        score = tot - score
        return score, None
