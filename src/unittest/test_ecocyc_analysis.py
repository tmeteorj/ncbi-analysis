import json
import os
import unittest
from urllib import request

from analysis.gene_name_eco_download import EcocycAnalysis
from utils.html_parser_util import EcocycHTMLParser, UrlHTMLParser


class TestEcocycAnalysis(unittest.TestCase):
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
        body = self.get_body('EG10875')
        parser = EcocycHTMLParser()
        parser.feed(body)
        self.assertEqual('cytosol, ribosome', parser.extract_attr['location'])
        self.assertEqual('rplN', parser.extract_attr['gene'])
        self.assertEqual('50S ribosomal subunit protein L14', parser.extract_attr['polypeptide'])
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

    def test_urls_extract(self):
        input_path = os.path.join(self.data_directory, 'transcription_ul.txt')
        with open(input_path, 'r', encoding='utf8') as fr:
            body = ''.join(fr.readlines())
        parser = UrlHTMLParser()
        parser.feed(body)
        for url, mock_name, title in parser.ecocycs:
            if mock_name is None:
                print(url)
        return parser.ecocycs

    def test_promoter_parser(self):
        promoter = '<b>Promoter:</b> ampCp (Experim. ev.)  <BR><b>Tr.Activators:</b> (BolA) <BR><b>Tr.Inhibitors:</b> (BolA) <BR><b>Tr.Start site:</b> 4,379,003'
        result = EcocycAnalysis.extract_attr(promoter)
        self.assertEqual('ampCp:4,379,003')
