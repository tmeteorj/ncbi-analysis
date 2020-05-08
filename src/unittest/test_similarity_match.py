import unittest

from analysis.gene_similarity_match import count_similarity


class TestSimilarityMatch(unittest.TestCase):
    def test_count_similarity(self):
        text_distance_similarity = count_similarity('text_distance', 100, 'ACGTACG', 'ACGACGT', 0)
        self.assertEqual(71, text_distance_similarity)
        offset_test = count_similarity('text_distance', 100, 'ACGTACG', 'GCTACGACGT', 3)
        self.assertEqual(text_distance_similarity, offset_test)

        match_similarity = count_similarity('char_match', 100, 'ACGTACG', 'ACGACGT', 0)
        self.assertEqual(42, match_similarity)
