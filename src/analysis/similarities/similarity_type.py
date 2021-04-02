from enum import Enum
from typing import Union


class SimilarityType(Enum):
    TextEdit = 0
    Direct = 1
    Consistency = 2
    Pattern = 3
    Blat = 4

    @staticmethod
    def from_object(similarity_type: Union[str, 'SimilarityType']):
        if isinstance(similarity_type, SimilarityType):
            return similarity_type
        type_str = similarity_type.lower()
        if type_str in ['text', 'textedit', 'textdistance', 'text_edit', 'text_distance']:
            return SimilarityType.TextEdit
        elif type_str in ['direct']:
            return SimilarityType.Direct
        elif type_str in ['consistency', 'consistent']:
            return SimilarityType.Consistency
        elif type_str in ['pattern', 'pat']:
            return SimilarityType.Pattern
        elif type_str in ['blat']:
            return SimilarityType.Blat
        else:
            raise ValueError(f"{similarity_type} not found in SimilarityType")
