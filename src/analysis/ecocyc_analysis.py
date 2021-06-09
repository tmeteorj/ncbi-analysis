import gzip
import json
import os
import traceback
from dataclasses import dataclass
from urllib import request

from utils.factories.logger_factory import LoggerFactory
from utils.gene_promoter_util import GeneTUInfo, get_target_promoter, get_all_promoters
from utils.html_parser_util import EcocycHTMLParser, UrlHTMLParser, GoHTMLParser
from utils.str_util import StrConverter


@dataclass
class EcocycAnalysis:
    input_path: str
    download_directory: str
    output_directory: str
    ecocyc_params: dict
    cookie: str = None

    def __post_init__(self):
        self.from_gene_names = self.ecocyc_params['from_gene_names']
        self.output_best_promoter = self.ecocyc_params['output_best_promoter']
        self.output_detail_information = self.ecocyc_params['output_detail_information']
        self.analysis_promoter = self.ecocyc_params['analysis_promoter']
        self.if_get_summary = self.ecocyc_params['if_get_summary']
        self.if_get_go_table = self.ecocyc_params['if_get_go_table']
        self.sequence_start_idx = None
        self.sequence_end_idx = None
        self.headers = {}
        self.inv_headers = []

        file_name = os.path.basename(self.input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.ecocyc_result_path = os.path.join(self.output_directory, '%s_ecocyc_result.txt' % file_prefix)
        self.ecocyc_error_path = os.path.join(self.output_directory, '%s_ecocyc_error.txt' % file_prefix)
        self.logger = LoggerFactory()

    def run(self):
        if self.from_gene_names:
            self.work_from_gene_list_file()
        else:
            self.work_from_url_list_file()

    def generate_header(self, items):
        for idx, col_name in enumerate(items.strip().split('\t')):
            self.headers[col_name] = idx
            self.inv_headers.append(col_name)
        self.sequence_end_idx = self.headers.get('gene_start_pos')
        self.sequence_start_idx = self.headers.get('promoter_pos')

    def work_from_gene_list_file(self):
        solve_cnt = 0
        succ_cnt = 0
        fail_cnt = 0
        fail_json_cnt = 0
        fw_error = open(self.ecocyc_error_path, 'w', encoding='utf8')
        fw_result = open(self.ecocyc_result_path, 'w', encoding='utf8')
        gene_items = list(
            filter(lambda arg: arg.strip() != '', open(self.input_path, 'r', encoding='utf8').readlines()))
        fw_result.write(gene_items[0])
        self.generate_header(gene_items[0])
        total_cnt = len(gene_items) - 1
        self.logger.info_with_expire_time(
            'Ecocyc analysis %d/%d=%.2f%%' % (solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt),
            solve_cnt, total_cnt)
        for line in gene_items[1:]:
            try:
                result = {}
                infos = line.strip().split('\t')
                for idx, info in enumerate(infos):
                    result[self.inv_headers[idx]] = info
                gene_name = result['gene']
                if gene_name.find('->') > 0:
                    gene_name, result['gene'] = result['gene'].split('->')
                self.write_body(gene_name=gene_name)
                ecocyc_id = self.get_ecocyc_id(prefix='gene_', gene_name=gene_name)
                result['ecocyc_id'] = ecocyc_id
                self.write_body(ecocyc_id=ecocyc_id, page_type="tu")
                self.analysis_xml(prefix='tu_', ecocyc_id=ecocyc_id, result=result)
                if self.if_get_summary:
                    self.write_body(ecocyc_id=ecocyc_id, page_type="summary")
                    self.analysis_xml(prefix='summary_', ecocyc_id=ecocyc_id, result=result)
                if self.analysis_promoter:
                    flag_json = self.write_body(ecocyc_id=ecocyc_id, page_type="promoter")
                    if flag_json:
                        self.analysis_json(prefix='promoter_', ecocyc_id=ecocyc_id, result=result,
                                           gene_name=result['gene'])
                    if not flag_json:
                        fail_json_cnt += 1
                if self.if_get_go_table:
                    self.write_body(ecocyc_id=ecocyc_id, page_type='go')
                    self.analysis_xml(prefix='go_', ecocyc_id=ecocyc_id, result=result)
                if result['gene'] != gene_name:
                    result['gene'] = gene_name + '->' + result['gene']
                fw_result.write(self.extract_output(result) + '\n')
                fw_result.flush()
                succ_cnt += 1
            except:
                fw_result.write('%s\tNot Found\n' % result['gene'])
                traceback.print_exc()
                fw_error.write(gene_name + '\n')
                fw_error.flush()
                fail_cnt += 1
            solve_cnt += 1
            self.logger.info_with_expire_time(
                'Ecocyc analysis %d/%d=%.2f%%, success_cnt=%d, fail_cnt=%d, json_download_fail=%d' % (
                    solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt,
                    succ_cnt, fail_cnt, fail_json_cnt),
                solve_cnt, total_cnt)
        fw_error.close()
        fw_result.close()

    def extract_output(self, result):
        output = []
        for name in self.inv_headers:
            if name == 'product_type':
                for key in ['enzyme', 'rna', 'protein', 'polypeptide', 'function when intact', 'transporter']:
                    if result.get(key, '') != '':
                        result['product_type'] = key
                        result['product'] = result[key]
            elif result.get(name, '') in ['', 'Not Found']:
                try:
                    if name in ['status', 'promoter_name', 'promoter_pos', 'gene_start_pos']:
                        if result['table_unites'][0] == 'Not Found':
                            if name == 'status':
                                result['status'] = 'Not Found'
                        else:
                            promoter = result['table_unites'][1]
                            result['status'] = 'Found'
                            result['gene_start_pos'] = result['table_unites'][0]
                            result['promoter_name'] = promoter.get_promoter_name()
                            result['promoter_pos'] = promoter.get_promoter_start_site(int_pos=True)
                except:
                    pass
            output.append(str(result.get(name, '')))
        return '\t'.join(output)

    def work_from_url_list_file(self):
        solve_cnt = 0
        succ_cnt = 0
        fail_cnt = 0
        fail_json_cnt = 0
        buff = []
        fw_error = open(self.ecocyc_error_path, 'w', encoding='utf8')
        fw_result = open(self.ecocyc_result_path, 'w', encoding='utf8')
        items = self.extract_urls_from_file()
        total_cnt = len(items)
        self.logger.info_with_expire_time(
            'Ecocyc analysis %d/%d=%.2f%%' % (solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt),
            solve_cnt, total_cnt)
        for url, mock_name, title in items:
            try:
                result = {}
                self.write_body(url=url, mock_name=mock_name)
                ecocyc_id = self.analysis_xml(prefix='url_', gene_name=mock_name, result=result)
                flag_json = False
                if ecocyc_id is not None:
                    result['ecocyc_id'] = ecocyc_id
                    flag_json = self.write_body(ecocyc_id=ecocyc_id, page_type="promoter")
                    if flag_json:
                        self.analysis_json(prefix='promoter_', ecocyc_id=ecocyc_id, result=result,
                                           gene_name=result['gene'])
                if not flag_json:
                    fail_json_cnt += 1
                temp = self.format_result_json(result, fw_error)
                if temp.strip() == '':
                    raise ValueError('No result found')
                buff.append(temp)
                fw_result.write(buff[-1])
                max_col = max(max_col, len(buff[-1].split('\t')))
                fw_result.flush()
                succ_cnt += 1
            except:
                traceback.print_exc()
                fw_error.write(url + '\t' + mock_name + '\t' + title + '\n')
                fw_error.flush()
                fail_cnt += 1
            solve_cnt += 1
            self.logger.info_with_expire_time(
                'Ecocyc analysis %d/%d=%.2f%%, success_cnt=%d, fail_cnt=%d, json_download_fail=%d' % (
                    solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt,
                    succ_cnt, fail_cnt, fail_json_cnt),
                solve_cnt, total_cnt)
        fw_error.close()
        fw_result.close()

    def extract_urls_from_file(self):
        with open(self.input_path, 'r', encoding='utf8') as fr:
            body = ''.join(fr.readlines())
        parser = UrlHTMLParser()
        parser.feed(body)
        return parser.ecocycs

    def write_body(self, url=None, mock_name=None, ecocyc_id=None, gene_name=None, page_type="tu"):
        if url is not None:
            urls = [url]
            origin_path = os.path.join(self.download_directory, mock_name + '.html')
            file_path = os.path.join(self.download_directory, 'url_' + mock_name + '.html')
            self.transform_file(origin_path, file_path)
        elif gene_name is not None:
            urls = ['http://ecocyc.org/ECOLI/substring-search?type=GENE&object=%s&geneSearch=Gene+Search' % gene_name]
            origin_path = os.path.join(self.download_directory, gene_name + '.html')
            file_path = os.path.join(self.download_directory, 'gene_' + gene_name + '.html')
            self.transform_file(origin_path, file_path)
        elif ecocyc_id is not None:
            if page_type == "tu":
                urls = ['https://ecocyc.org/gene?orgid=ECOLI&id=%s#tab=TU' % ecocyc_id]
                origin_path = os.path.join(self.download_directory, ecocyc_id + '.html')
                file_path = os.path.join(self.download_directory, 'tu_' + ecocyc_id + '.html')
                self.transform_file(origin_path, file_path)
            elif page_type == "promoter":
                urls = [

                    'https://ecocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://ecocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://ecocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-INDEX_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://ecocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-INDEX_NO-PLOC_%s.wg' % ecocyc_id
                ]
                origin_path = os.path.join(self.download_directory, ecocyc_id + '.json')
                file_path = os.path.join(self.download_directory, 'promoter_' + ecocyc_id + '.json')
                self.transform_file(origin_path, file_path)
            elif page_type == "summary":
                urls = ['https://biocyc.org/gene-tab?id=%s&orgid=ECOLI&tab=SUMMARY' % ecocyc_id]
                file_path = os.path.join(self.download_directory, 'summary_' + ecocyc_id + '.html')
            elif page_type == "go":
                urls = ['https://biocyc.org/gene-tab?id=%s&orgid=ECOLI&tab=GO' % ecocyc_id]
                file_path = os.path.join(self.download_directory, 'go_' + ecocyc_id + '.html')
        else:
            raise ValueError('Parameter not correct')
        if os.path.exists(file_path):
            return True
        headers = {"Host": "ecocyc.org",
                   "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.108 Safari/537.36",
                   "Accept": "*/*",
                   "Sec-Fetch-Site": "same-origin",
                   "Sec-Fetch-Mode": "cors",
                   'Accept-Encoding': "gzip, deflate, br",
                   'Connection': "Keep-Alive",
                   'Cookie': self.cookie
                   }
        for retry_time in range(3):
            flag = False
            for url in urls:
                try:
                    if retry_time == 0:
                        x = request.urlopen(url, timeout=30)
                        body = x.read().decode('utf8')
                        with open(file_path, 'w', encoding='utf8') as fw:
                            fw.write(body)
                            flag = True
                            break
                    elif retry_time == 1:
                        url = "https://biocyc.org/tmp/ptools-images/ECOLI/%s_REG-SUMMARY.wg" % ecocyc_id
                        req = request.Request(url=url, headers=headers)
                        x = request.urlopen(req, timeout=30)
                        body = x.read()
                        break
                    else:
                        req = request.Request(url=url, headers=headers)
                        x = request.urlopen(req, timeout=30)
                        body = x.read()
                        for item in x.headers._headers:
                            if item[0].lower() == 'content-encoding' and item[1].lower() == 'gzip':
                                body = gzip.decompress(body)
                        body = body.decode('utf-8')
                        with open(file_path, 'w', encoding='utf8') as fw:
                            fw.write(body)
                            flag = True
                            break
                except:
                    continue
            if flag:
                break
        return flag

    def analysis_xml(self, prefix, ecocyc_id, result):
        xml_path = os.path.join(self.download_directory, prefix + ecocyc_id + '.html')
        with open(xml_path, 'r', encoding='utf8') as fr:
            body = ''.join(fr.readlines())
        if prefix == 'summary_':
            parser = EcocycHTMLParser(do_extract_summary=True)
            parser.feed(''.join(body))
            result['summary'] = parser.extract_attr['summary']
        elif prefix == 'go_':
            parser = GoHTMLParser()
            parser.feed(''.join(body))
            result['go'] = ';'.join(['%s=%s' % (k, v) for k, v in parser.go_table])
        else:
            parser = EcocycHTMLParser()
            parser.feed(''.join(body))
            for k, v in parser.extract_attr.items():
                if k == 'map position':
                    result['map_start_pos'] = v[0]
                    result['map_end_pos'] = v[1]
                elif v is not None:
                    result[k] = v.strip('__#####__')
            return parser.ecocyc_id

    def analysis_json(self, prefix, ecocyc_id, result, gene_name=None):
        json_path = os.path.join(self.download_directory, prefix + ecocyc_id + '.json')
        with open(json_path, 'r') as fr:
            body = ''.join(fr.readlines())
        body = json.loads(body)
        data = []
        target_gene = None
        for link in body['links']:
            gene_tu = GeneTUInfo(link)
            if self.output_best_promoter and gene_name is not None:
                if gene_tu.is_gene(gene_name):
                    target_gene = gene_tu
            data.append(gene_tu)
        if self.output_best_promoter:
            flag = False
            if target_gene is not None:
                target_promoter, near_gene_pos = get_target_promoter(target_gene, data)
                if target_promoter is not None:
                    data = [near_gene_pos, target_promoter]
                    flag = True
            if not flag:
                data = ['Not Found']
        else:
            data = get_all_promoters(data, True)
        result['table_unites'] = data

    def get_ecocyc_id(self, prefix, gene_name):
        xml_path = os.path.join(self.download_directory, prefix + gene_name + '.html')
        with open(xml_path, 'r') as fr:
            body = ''.join(fr.readlines())

        parser = EcocycHTMLParser(do_extract_id=True, gene_name=gene_name)
        parser.feed(''.join(body))
        if parser.ecocyc_id is None:
            raise RuntimeError('Ecocyc is is None, parse error for %s' % gene_name)
        return parser.ecocyc_id

    @staticmethod
    def transform_file(original_path, new_path):
        if not os.path.exists(new_path) and os.path.exists(original_path):
            os.rename(original_path, new_path)

    def format_result_json(self, result, fw_error=None):
        keys = ['gene', 'cluster']
        info = [result.get(key, '') for key in keys]
        product_type = ''
        product = ''
        for key in ['rna', 'protein', 'polypeptide', 'enzyme', 'function when intact', 'transporter']:
            val = result.get(key, '')
            if val is None or val == '': continue
            product_type = key
            product = val
        info.extend([product_type, product])
        for key in ['location']:  # , 'reaction']:
            val = result.get(key, '')
            if val is None: val = ''
            info.append(val)
        table_unites = result.get('table_unites', [])
        if self.output_best_promoter and len(table_unites) == 2 and type(table_unites[0]) is int:
            near_gene_pos, promoter = table_unites
            info.extend(['Gene Start Position', near_gene_pos])
            info.extend([promoter.get_promoter_name(), promoter.get_promoter_start_site(int_pos=True)])
            info = list(map(str, info))
        else:
            for promoter in table_unites:
                info.extend([promoter.get_promoter_name(), promoter.get_promoter_start_site()])
            info = list(map(str, info))
            if self.output_best_promoter and fw_error is not None:
                fw_error.write('\t'.join(info) + '\n')
        return '\t'.join(info) + '\n'
