from dataclasses import field, dataclass
from typing import Mapping

from analysis.models.similarity_type import SimilarityType


@dataclass
class MatchCandidate:
    left: int
    right: int
    is_reverse: bool
    database_length: int
    weighted_similarity: float
    similarity_dict: Mapping[SimilarityType, float] = field(default_factory=dict)

    def __post_init__(self):
        if self.is_reverse:
            self.start = self.database_length - self.left
            self.end = self.database_length - self.right
        else:
            self.start = self.left + 1
            self.end = self.right + 1
        self.should_ignore = False
        self.original_match_left = self.left
        self.original_match_right = self.right

    def get_similarity_str(self):
        result = 'weighted=%.2f' % self.weighted_similarity
        for key, value in self.similarity_dict.items():
            result += ', %s=%.2f' % (key.name, value)
        return result

    def __str__(self):
        return '[%d-%d], (%s), %s, ' % (self.start, self.right, self.get_similarity_str(), self.should_ignore)

    def __le__(self, other):
        return self.weighted_similarity <= other.weighted_similarity

    def __lt__(self, other):
        return self.weighted_similarity < other.weighted_similarity

    def __ge__(self, other):
        return self.weighted_similarity >= other.weighted_similarity

    def __gt__(self, other):
        return self.weighted_similarity > other.weighted_similarity
