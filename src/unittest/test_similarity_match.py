import unittest

from analysis.gene_similarity_match import count_similarity, fast_skip, count_acgt, compute_blat_similarity


class TestSimilarityMatch(unittest.TestCase):
    def test_count_similarity(self):
        text_distance_similarity = count_similarity('text_distance', 100, 'ACGTACG', 'ACGACGT', 0, is_reverse=False)
        self.assertEqual(71, text_distance_similarity)
        offset_test = count_similarity('text_distance', 100, 'ACGTACG', 'GCTACGACGT', 3, is_reverse=False)
        self.assertEqual(text_distance_similarity, offset_test)

        match_similarity = count_similarity('char_match', 100, 'ACGTACG', 'ACGACGT', 0, is_reverse=False)
        self.assertEqual(42, match_similarity)

    def test_fast_skip(self):
        source_gene = 'AAAATTTAA'
        target_gene = 'AAATTTTGG'
        pat = None
        source_count = count_acgt(source_gene)
        self.assertFalse(fast_skip(source_count, 9, target_gene, 0, 3, pat))
        self.assertTrue(fast_skip(source_count, 9, target_gene, 0, 7, pat))
        pat = '.*AA.*GG.*'
        self.assertFalse(fast_skip(source_count, 9, target_gene, 0, 3, pat))
        pat = '.*AA.*AA'
        self.assertTrue(fast_skip(source_count, 9, target_gene, 0, 3, pat))

    def test_blat_similarity(self):
        gene = "tgatatca"
        test_cases = [
            ["atgatatca", 0, False, None],
            ["tgatatca", 0, False, None],
            ["tgataatca", 0, True, 9],
            ["tagataatca", 0, True, 10],
            ["tagtactaatca", 0, False, None],
            ["taaaagataatca", 0, True, 13]
        ]
        for database, offset, expect_flag, expect_pos in test_cases:
            flag, pos = compute_blat_similarity(gene, database, offset)
            self.assertEqual(expect_flag, flag)
            self.assertEqual(expect_pos, pos)
