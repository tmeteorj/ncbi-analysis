from typing import Union

from analysis.similarities.blat_similarity import BlatSimilarity
from analysis.similarities.consistency_similarity import ConsistencySimilarity
from analysis.similarities.direct_match_similarity import DirectMatchSimilarity
from analysis.similarities.pattern_similarity import PatternSimilarity, MatchPattern
from analysis.enum_types import SimilarityType
from analysis.similarities.text_edit_similarity import TextEditSimilarity


class SimilarityFactory(object):
    @staticmethod
    def get_similarity(
            similarity_type: Union[str, SimilarityType],
            continuous_mismatch_limit: int = 10,
            max_patience: int = 2,
            match_pattern: MatchPattern = None,
            mid_limit: int = 10,
            end_limit: int = 2
    ):
        similarity_type = SimilarityType.from_object(similarity_type)
        if similarity_type == SimilarityType.TextEdit:
            return TextEditSimilarity(continuous_mismatch_limit=continuous_mismatch_limit)
        elif similarity_type == SimilarityType.Direct:
            return DirectMatchSimilarity()
        elif similarity_type == SimilarityType.Consistency:
            return ConsistencySimilarity(max_patience=max_patience)
        elif similarity_type == SimilarityType.Pattern:
            return PatternSimilarity(match_pattern=match_pattern)
        elif similarity_type == SimilarityType.Blat:
            return BlatSimilarity(mid_limit=mid_limit, end_limit=end_limit)
        else:
            raise ValueError(f"Not supported type {similarity_type}")
