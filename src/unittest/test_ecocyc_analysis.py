import json
import os
import unittest
from urllib import request

from utils.gene_promoter_util import GeneTUInfo, filter_same_direction, get_all_promoters, get_target_promoter
from utils.html_parser_util import EcocycHTMLParser, UrlHTMLParser


class TestEcocycAnalysis(unittest.TestCase):
    def setUp(self):
        self.root_directory = os.sep.join(os.getcwd().split(os.sep)[:-2])
        self.data_directory = os.path.join(self.root_directory, 'data', 'rna_analysis')
        self.download_directory = os.path.join(self.data_directory, 'ecocyc_download_data')

    def get_body(self, gene_name, prefix, suffix='.html'):
        path = os.path.join(self.download_directory, prefix + gene_name + suffix)
        with open(path, 'r') as f:
            body = ''.join(f.readlines())
        return body

    def test_simple_ecocyc_id_extract(self):
        body = self.get_body('fdhE', prefix='gene_')
        parser = EcocycHTMLParser(do_extract_id=True)
        parser.feed(body)
        self.assertEqual('EG10284', parser.ecocyc_id)

    def test_multi_ecocyc_id_extract(self):
        body = self.get_body('fadI', prefix='gene_')
        parser = EcocycHTMLParser(do_extract_id=True, gene_name='fadI')
        parser.feed(body)
        self.assertEqual('G7213', parser.ecocyc_id)

    def test_attr_extraction(self):
        body = self.get_body('EG10875', 'tu_')
        parser = EcocycHTMLParser()
        parser.feed(body)
        self.assertEqual('cytosol, ribosome', parser.extract_attr['location'])
        self.assertEqual('rplN', parser.extract_attr['gene'])
        self.assertEqual('50S ribosomal subunit protein L14', parser.extract_attr['polypeptide'])
        # self.assertIsNotNone(parser.extract_attr['Evidence'])

    def test_json_download(self):
        ecocyc_id = 'EG10285'
        url = "https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-INDEX_NO-PLOC_%s.wg" % ecocyc_id
        file_path = os.path.join(self.download_directory, 'promoter_' + ecocyc_id + '.json')
        if not os.path.exists(file_path):
            x = request.urlopen(url, timeout=60)
            body = x.read().decode('utf8')
            with open(file_path, 'w', encoding='utf8') as fw:
                fw.write(body)
        body = json.loads(self.get_body(ecocyc_id, 'promoter_', '.json'))
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

    def test_promoter_position_extract(self):
        gene_name = 'napF'
        input_path = os.path.join(self.download_directory, 'promoter_EG12068.json')
        json_path = os.path.join(input_path)
        with open(json_path, 'r') as fr:
            body = ''.join(fr.readlines())
        body = json.loads(body)
        data = []
        gene_tu = None
        for link in body['links']:
            gene_tu_info = GeneTUInfo(link)
            if gene_tu_info.is_gene(gene_name):
                gene_tu = gene_tu_info
            data.append(gene_tu_info)
        self.assertEqual(39, len(data))
        self.assertEqual(15, gene_tu.idx)
        promoters = get_all_promoters(data)
        self.assertEqual(7, len(promoters))
        same_direction_promoters = filter_same_direction(gene_tu, promoters)
        self.assertEqual(5, len(same_direction_promoters))
        promoter, _ = get_target_promoter(gene_tu, data)
        self.assertEqual(10, promoter.idx)

    def test_promoter_position_extract_with_other_gene(self):
        gene_name = 'tpiA'
        input_path = os.path.join(self.download_directory, 'promoter_EG11015.json')
        json_path = os.path.join(input_path)
        with open(json_path, 'r') as fr:
            body = ''.join(fr.readlines())
        body = json.loads(body)
        data = []
        gene_tu = None
        for link in body['links']:
            gene_tu_info = GeneTUInfo(link)
            if gene_tu_info.is_gene(gene_name):
                gene_tu = gene_tu_info
            data.append(gene_tu_info)
        self.assertEqual(10, len(data))
        self.assertEqual(1, gene_tu.idx)
        promoters = get_all_promoters(data)
        self.assertEqual(3, len(promoters))
        same_direction_promoters = filter_same_direction(gene_tu, promoters)
        self.assertEqual(3, len(same_direction_promoters))
        promoter, _ = get_target_promoter(gene_tu, data)
        self.assertEqual(0, promoter.idx)

    def test_target_promoter(self):
        test_cases = [['rplJ', 'EG10871', 0],
                      ['rplA', 'EG10864', 1],
                      ['ftnB', 'G7033', 2]
                      ]
        for gene_name, ecocyc_id, target_idx in test_cases:
            body = self.get_body(ecocyc_id, 'promoter_', '.json')
            body = json.loads(body)
            data = []
            target_gene = None
            for link in body['links']:
                gene_tu_info = GeneTUInfo(link)
                if gene_tu_info.is_gene(gene_name):
                    target_gene = gene_tu_info
                data.append(gene_tu_info)
            promoter, _ = get_target_promoter(target_gene, data)
            self.assertEqual(target_idx, promoter.idx)

    def test_target_promoter_with_position(self):
        test_cases = [['nuoL', 'EG12092', 'nuoAp1', 2405072],
                      ['bipA', 'EG11837', 'typAp1', 4058407]]
        for gene_name, ecocyc_id, target_promoter_name, target_pos in test_cases:
            body = self.get_body(ecocyc_id, 'promoter_', '.json')
            body = json.loads(body)
            data = []
            target_gene = None
            for link in body['links']:
                gene_tu_info = GeneTUInfo(link)
                if gene_tu_info.is_gene(gene_name):
                    target_gene = gene_tu_info
                data.append(gene_tu_info)
            promoter, near_gene_pos = get_target_promoter(target_gene, data)
            self.assertTrue(promoter.get_promoter_name().startswith(target_promoter_name),
                            'expect= %s, actual= %s' % (target_promoter_name, promoter.get_promoter_name()))
            self.assertEqual(target_pos, near_gene_pos)
