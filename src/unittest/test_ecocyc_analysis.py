import gzip
import json
import os
import unittest
from urllib import request

from analysis.ecocyc_analysis import EcocycAnalysis
from utils.gene_promoter_util import GeneTUInfo, filter_same_direction, get_all_promoters, get_target_promoter
from utils.html_parser_util import EcocycHTMLParser, UrlHTMLParser


class TestEcocycAnalysis(unittest.TestCase):
    def setUp(self):
        self.root_directory = os.sep.join(os.getcwd().split(os.sep)[:-2])
        self.data_directory = os.path.join(self.root_directory, 'data', 'rna_analysis')
        self.download_directory = os.path.join(self.data_directory, 'ecocyc_download_data')
        self.ecocyc_analysis = EcocycAnalysis('',
                                              self.download_directory,
                                              os.path.join(self.root_directory, 'data', 'rna_analysis_result'),
                                              {
                                                  'from_gene_names': True,
                                                  'output_best_promoter': True,
                                                  'output_gene_sequence': True,
                                                  'output_detail_information': True,
                                                  'analysis_promoter': True,
                                                  'if_get_summary': True
                                              })

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
        for name, ecocyc_id in [['fadI', 'G7213'],
                                ['sra', 'EG11508']]:
            body = self.get_body(name, prefix='gene_')
            parser = EcocycHTMLParser(do_extract_id=True, gene_name=name)
            parser.feed(body)
            self.assertEqual(ecocyc_id, parser.ecocyc_id)

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

    def test_target_promoter_with_gene(self):
        test_cases = [['rpsR', 'rpsFp', 4425118],
                      ['nfsB', 'nfsBp', 605424]]
        for target_gene, target_promoter_name, target_pos in test_cases:
            result = {}
            self.ecocyc_analysis.write_body(gene_name=target_gene)
            ecocyc_id = self.ecocyc_analysis.get_ecocyc_id(prefix='gene_', gene_name=target_gene)
            self.assertIsNotNone(ecocyc_id)
            self.ecocyc_analysis.write_body(ecocyc_id=ecocyc_id, page_type=True)
            flag_json = self.ecocyc_analysis.write_body(ecocyc_id=ecocyc_id, page_type=False)
            self.assertTrue(flag_json)
            _ = self.ecocyc_analysis.analysis_xml(prefix='tu_', ecocyc_id=ecocyc_id, result=result)
            self.ecocyc_analysis.analysis_json(prefix='promoter_', ecocyc_id=ecocyc_id, result=result,
                                               gene_name=result['gene'])
            table_unites = result['table_unites']
            self.assertEqual(2, len(table_unites))
            self.assertEqual(target_pos, table_unites[0])
            self.assertTrue(table_unites[1].get_promoter_name().startswith(target_promoter_name), target_gene)

    def test_json_download(self):
        url = 'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-INDEX_NO-PLOC_G6081.wg'
        headers = {"Host": "biocyc.org",
                   "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.108 Safari/537.36",
                   "Accept": "*/*",
                   "Sec-Fetch-Site": "same-origin",
                   "Sec-Fetch-Mode": "cors",
                   "Referer": "https://biocyc.org/gene?orgid=ECOLI&id=EG10917",
                   'Accept-Encoding': "gzip, deflate, br",
                   'Connection': "Keep-Alive",
                   'Cookie': 'pagecount=6; frameHeight=761; frameWidth=1500; PTools-session=biocyc14b~biocyc14-3786098971%7CNIL%20NIL%20%22%22%20NIL%200%20(%3AWEB%20NIL%20-1%20((%3ABASICS%20-1)%20(%3AQUERIES%20-1)%20(%3AADVANCED%20-1)))%20NIL%20NIL%20ECOBASE%20NIL%20NIL%20%7Ch0x9lugmbx6fhk3bevsf98cen2omar5; _ga=GA1.2.407871027.1577110083; JSESSIONID=C648FCE2DAC5362DE3EE9322D39BD47E'
                   }
        req = request.Request(url=url, headers=headers)
        x = request.urlopen(req, timeout=30)
        body = x.read()
        dec_body = gzip.decompress(body)
        dec_body = dec_body.decode('utf-8')
        print(body)

    def test_summary_extract(self):
        body = self.get_body('EG30113', 'summary_')
        parser = EcocycHTMLParser(do_extract_summary=True)
        parser.feed(body)
        summary = parser.extract_attr['summary']
        self.assertIsNotNone(summary)
