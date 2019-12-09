import json
import os
import traceback
from urllib import request

from utils.factories.logger_factory import LoggerFactory
from utils.gene_promoter_util import GeneTUInfo, get_target_promoter, get_all_promoters
from utils.html_parser_util import EcocycHTMLParser, UrlHTMLParser
from utils.str_util import StrConverter


class EcocycAnalysis:
    def __init__(self, input_path, download_directory, output_directory, from_gene_names=True,
                 output_best_promoter=False):
        self.download_directory = download_directory
        self.output_directory = output_directory
        self.input_path = input_path
        self.from_gene_names = from_gene_names
        self.output_best_promoter = output_best_promoter

        file_name = os.path.basename(input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.ecocyc_result_path = os.path.join(self.output_directory, '%s_ecocyc_result.txt' % file_prefix)
        self.ecocyc_error_path = os.path.join(self.output_directory, '%s_ecocyc_error.txt' % file_prefix)
        self.logger = LoggerFactory()

    def run(self):
        solve_cnt = 0
        succ_cnt = 0
        fail_cnt = 0
        fail_json_cnt = 0
        fw_error = open(self.ecocyc_error_path, 'w', encoding='utf8')
        fw_result = open(self.ecocyc_result_path, 'w', encoding='utf8')
        buff = []
        max_col = 0
        if self.from_gene_names:
            gene_items = list(filter(lambda arg: arg.strip() != '', open(self.input_path, 'r').readlines()))
            total_cnt = len(gene_items)
            self.logger.info_with_expire_time(
                'Ecocyc analysis %d/%d=%.2f%%' % (solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt),
                solve_cnt, total_cnt)
            for line in gene_items:
                try:
                    info = line.strip().split('\t')
                    gene_name = info[0].strip()
                    result = {'gene': gene_name}
                    if len(info) > 1:
                        result['cluster'] = info[1]
                    self.write_body(gene_name=gene_name)
                    ecocyc_id = self.get_ecocyc_id(prefix='gene_', gene_name=gene_name)
                    result['ecocyc_id'] = ecocyc_id
                    self.write_body(ecocyc_id=ecocyc_id, get_summary=True)
                    flag_json = self.write_body(ecocyc_id=ecocyc_id, get_summary=False)
                    _ = self.analysis_xml(prefix='tu_', ecocyc_id=ecocyc_id, result=result)
                    if flag_json:
                        self.analysis_json(prefix='promoter_', ecocyc_id=ecocyc_id, result=result,
                                           gene_name=result['gene'])
                    if not flag_json:
                        fail_json_cnt += 1
                    if result['gene'] != gene_name:
                        result['gene'] = gene_name + '->' + result['gene']
                    buff.append(self.format_result_json(result, fw_error))
                    fw_result.write(buff[-1])
                    max_col = max(max_col, len(buff[-1].split('\t')))
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
        else:
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
                        flag_json = self.write_body(ecocyc_id=ecocyc_id, get_summary=False)
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
        with open(self.ecocyc_result_path, 'w', encoding='utf8') as fw:
            fw.write('gene\tcluster\tproduct_type\tproduct\tlocation\treaction')
            max_col -= 6
            max_col //= 2
            for idx in range(max_col):
                fw.write('\tpromoter#%d\tstart#%d' % (idx + 1, idx + 1))
            fw.write('\n')
            for line in buff:
                if line.strip() == '': continue
                fw.write(line)

    def extract_urls_from_file(self):
        with open(self.input_path, 'r', encoding='utf8') as fr:
            body = ''.join(fr.readlines())
        parser = UrlHTMLParser()
        parser.feed(body)
        return parser.ecocycs

    def write_body(self, url=None, mock_name=None, ecocyc_id=None, gene_name=None, get_summary=True):
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
            if get_summary:
                urls = ['https://ecocyc.org/gene?orgid=ECOLI&id=%s#tab=TU' % ecocyc_id]
                origin_path = os.path.join(self.download_directory, ecocyc_id + '.html')
                file_path = os.path.join(self.download_directory, 'tu_' + ecocyc_id + '.html')
                self.transform_file(origin_path, file_path)
            else:
                urls = [
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=-1_NO-INDEX_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-INDEX_NO-PLOC_%s.wg' % ecocyc_id,
                    'https://biocyc.org/tmp/ptools-images/ECOLI/TU_dir=1_topdir=1_NO-PLOC_%s.wg' % ecocyc_id
                ]
                origin_path = os.path.join(self.download_directory, ecocyc_id + '.json')
                file_path = os.path.join(self.download_directory, 'promoter_' + ecocyc_id + '.json')
                self.transform_file(origin_path, file_path)
        else:
            raise ValueError('Parameter not correct')
        if os.path.exists(file_path):
            return True
        for retry_time in range(1):
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

    def analysis_xml(self, prefix, ecocyc_id, result):
        xml_path = os.path.join(self.download_directory, prefix + ecocyc_id + '.html')
        with open(xml_path, 'r', encoding='utf8') as fr:
            body = ''.join(fr.readlines())
        parser = EcocycHTMLParser()
        parser.feed(''.join(body))
        for k, v in parser.extract_attr.items():
            if v is not None:
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
        if self.output_best_promoter and target_gene is not None:
            target_promoter = get_target_promoter(target_gene, data)
            if target_promoter is not None:
                data = [target_gene, target_promoter]
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
        if self.output_best_promoter and len(table_unites) == 2 and table_unites[0].is_gene():
            gene, promoter = table_unites
            info.extend(['Gene Start Position', gene.get_gene_start_position()])
            info.extend([promoter.get_promoter_name(), promoter.get_promoter_start_site(int_pos=True)])
            info = list(map(str, info))
        else:
            for promoter in table_unites:
                info.extend([promoter.get_promoter_name(), promoter.get_promoter_start_site()])
            info = list(map(str, info))
            if self.output_best_promoter and fw_error is not None:
                fw_error.write('\t'.join(info) + '\n')
        return '\t'.join(info) + '\n'
