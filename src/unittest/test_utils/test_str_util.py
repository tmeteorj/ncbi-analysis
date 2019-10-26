import unittest

from utils.str_util import StrConverter


class TestStrUtil(unittest.TestCase):

    def test_extract_filename(self):
        check_list = [
            ('18 rna utr', '18_rna_utr'),
            ('18_rna_utr.txt', '18_rna_utr')
        ]
        for in_name, expect_name in check_list:
            out_name = StrConverter.extract_file_name(in_name)
            self.assertEqual(expect_name, out_name, 'Extract file name failed: %s' % in_name)
