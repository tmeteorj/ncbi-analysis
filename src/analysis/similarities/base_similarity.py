from typing import Tuple, Any


class BaseSimilarity(object):

    def get_similarity(self, gene: str, database: str, offset: int) -> Tuple[float, Any]:
        raise NotImplementedError()

    def rendering_sequence(self, gene: str, database: str, offset: int) -> Tuple[list, list, list]:
        raise NotImplementedError()

    @staticmethod
    def should_change(a, b):
        if a == b:
            return 0
        elif a == 'c' and b == 't':
            return 0
        return 1
