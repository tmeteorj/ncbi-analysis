import json
import os
import unittest
from urllib import request

from utils.html_parser_util import EcocycHTMLParser


class TestHTMLParser(unittest.TestCase):
    def setUp(self):
        self.root_directory = os.sep.join(os.getcwd().split(os.sep)[:-2])
        self.data_directory = os.path.join(self.root_directory, 'data', 'rna_analysis')
        self.download_directory = os.path.join(self.data_directory, 'ecocyc_download_data')

    def get_body(self, gene_name, suffix='.html'):
        path = os.path.join(self.download_directory, gene_name + suffix)
        with open(path, 'r') as f:
            body = ''.join(f.readlines())
        return body

    def test_simple_ecocyc_id_extract(self):
        body = self.get_body('fdhE')
        parser = EcocycHTMLParser(do_extract_id=True)
        parser.feed(body)
        self.assertEqual('EG10284', parser.ecocyc_id)

    def test_multi_ecocyc_id_extract(self):
        body = self.get_body('fadI')
        parser = EcocycHTMLParser(do_extract_id=True, gene_name='fadI')
        parser.feed(body)
        self.assertEqual('G7213', parser.ecocyc_id)

    def test_attr_extraction(self):
        body = self.get_body('EG11227')
        parser = EcocycHTMLParser()
        parser.feed(body)
        self.assertIsNotNone(parser.extract_attr['Location'])
        # self.assertIsNotNone(parser.extract_attr['Evidence'])

    def test_json_download(self):
        ecocyc_id = 'EG10285'
        url = "https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-INDEX_NO-PLOC_%s.wg" % ecocyc_id
        file_path = os.path.join(self.download_directory, ecocyc_id + '.json')
        if not os.path.exists(file_path):
            x = request.urlopen(url, timeout=60)
            body = x.read().decode('utf8')
            with open(file_path, 'w', encoding='utf8') as fw:
                fw.write(body)
        body = json.loads(self.get_body(ecocyc_id, '.json'))
        for link in body['links']:
            attr = link[6]
            if attr.find('fdhFp') > 0:
                print(attr)
