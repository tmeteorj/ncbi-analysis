import json
import os
import re
import traceback
from urllib import request

from utils.factories.logger_factory import LoggerFactory
from utils.html_parser_util import EcocycHTMLParser
from utils.str_util import StrConverter


class EcocycAnalysis:
    def __init__(self, input_path, download_directory, output_directory):
        self.download_directory = download_directory
        self.output_directory = output_directory
        self.input_path = input_path

        file_name = os.path.basename(input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.ecocyc_result_path = os.path.join(self.output_directory, '%s_ecocyc_result.txt' % file_prefix)
        self.ecocyc_error_path = os.path.join(self.output_directory, '%s_ecocyc_error.txt' % file_prefix)
        self.logger = LoggerFactory()

    def run(self):
        gene_names = list(filter(lambda arg: arg.strip() != '', open(self.input_path, 'r').readlines()))
        total_cnt = len(gene_names)
        solve_cnt = 0
        succ_cnt = 0
        fail_cnt = 0
        fail_json_cnt = 0
        self.logger.info_with_expire_time(
            'Ecocyc analysis %d/%d=%.2f%%' % (solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt),
            solve_cnt, total_cnt)
        fw_error = open(self.ecocyc_error_path, 'w', encoding='utf8')
        fw_result = open(self.ecocyc_result_path, 'w', encoding='utf8')
        for gene_name in gene_names:
            try:
                gene_name = gene_name.strip()
                flag_id = self.write_body(gene_name=gene_name)
                ecocyc_id = self.get_ecocyc_id(gene_name)
                flag_xml = self.write_body(ecocyc_id=ecocyc_id, get_summary=True)
                flag_json = self.write_body(ecocyc_id=ecocyc_id, get_summary=False)
                xml_data = self.analysis_xml(ecocyc_id)
                json_data = self.analysis_json(gene_name, ecocyc_id) if flag_json else ''
                if not flag_json:
                    fail_json_cnt += 1
                fw_result.write('%s\t%s\t%s\t%s\n' % (gene_name, ecocyc_id, xml_data, json_data))
                fw_result.flush()
                succ_cnt += 1
            except:
                traceback.print_exc()
                fw_error.write(gene_name + '\t' + str(ecocyc_id) + '\n')
                fw_error.flush()
                fail_cnt += 1
            solve_cnt += 1
            self.logger.info_with_expire_time(
                'Ecocyc analysis %d/%d=%.2f%%, success_cnt=%d, fail_cnt=%d, json_download_fail=%d' % (
                    solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt,
                    succ_cnt, fail_cnt, fail_json_cnt),
                solve_cnt, total_cnt)
        fw_error.close()

    def write_body(self, ecocyc_id=None, gene_name=None, get_summary=True):
        if gene_name is not None:
            urls = ['http://ecocyc.org/ECOLI/substring-search?type=GENE&object=%s&geneSearch=Gene+Search' % gene_name]
            file_path = os.path.join(self.download_directory, gene_name + '.html')
        elif ecocyc_id is not None:
            if get_summary:
                urls = ['https://ecocyc.org/gene?orgid=ECOLI&id=%s#tab=TU' % ecocyc_id]
                file_path = os.path.join(self.download_directory, ecocyc_id + '.html')
            else:
                urls = [
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-INDEX_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-INDEX_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-PLOC_%s.wg' % ecocyc_id
                ]
                file_path = os.path.join(self.download_directory, ecocyc_id + '.json')
        else:
            raise ValueError('Parameter not correct')
        if os.path.exists(file_path):
            return True
        for retry_time in range(3):
            flag = False
            failed_urls = []
            for url in urls:
                try:
                    x = request.urlopen(url, timeout=30)
                    body = x.read().decode('utf8')
                    with open(file_path, 'w', encoding='utf8') as fw:
                        fw.write(body)
                        flag = True
                    break
                except:
                    failed_urls.append(url)
                    continue
            if flag:
                break
        return flag

    def analysis_xml(self, ecocyc_id):
        xml_path = os.path.join(self.download_directory, ecocyc_id + '.html')
        with open(xml_path, 'r', encoding='utf8') as fr:
            body = ''.join(fr.readlines())
        parser = EcocycHTMLParser()
        parser.feed(''.join(body))
        location = parser.extract_attr['Location']
        evidence = parser.extract_attr['Reaction']
        if location is None: location = ''
        if evidence is None: evidence = ''
        return '%s\t%s' % (re.sub(r'\s+', ' ', location), re.sub(r'\s+', ' ', evidence))

    def analysis_json(self, gene_name, ecocyc_id):
        json_path = os.path.join(self.download_directory, ecocyc_id + '.json')
        with open(json_path, 'r') as fr:
            body = ''.join(fr.readlines())
        body = json.loads(body)
        data = []
        for link in body['links']:
            attr = link[-1]
            if attr.find('Tr.Start') > 0:
                data.append(self.extract_attr(attr))
        return '\t'.join(data)

    def extract_attr(self, attr):
        attr = attr.replace('<b>', '').replace('</b>', '')
        return re.sub(r'\s+', ' ', re.sub(r'\s+<BR>\s+', ';', attr))

    def get_ecocyc_id(self, gene_name):
        xml_path = os.path.join(self.download_directory, gene_name + '.html')
        with open(xml_path, 'r') as fr:
            body = ''.join(fr.readlines())

        parser = EcocycHTMLParser(do_extract_id=True, gene_name=gene_name)
        parser.feed(''.join(body))
        if parser.ecocyc_id is None:
            raise RuntimeError('Ecocyc is is None, parse error for %s' % gene_name)
        return parser.ecocyc_id
