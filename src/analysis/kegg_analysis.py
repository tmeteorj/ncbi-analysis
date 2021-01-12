import gzip
import json
import os
import traceback
from dataclasses import dataclass
from urllib import request

from utils.factories.logger_factory import LoggerFactory
from utils.gene_promoter_util import GeneTUInfo, get_target_promoter, get_all_promoters
from utils.html_parser_util import EcocycHTMLParser, UrlHTMLParser, GoHTMLParser, KeggIdHTMLParser, \
    KeggPathwayHTMLParser
from utils.str_util import StrConverter


@dataclass
class KeggAnalysis:
    input_path: str
    download_directory: str
    output_directory: str
    is_gene: bool = True

    def __post_init__(self):
        os.makedirs(self.download_directory, exist_ok=True)
        file_name = os.path.basename(self.input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.kegg_result_path = os.path.join(self.output_directory, '%s_kegg_result.txt' % file_prefix)
        self.kegg_error_path = os.path.join(self.output_directory, '%s_kegg_error.txt' % file_prefix)
        self.logger = LoggerFactory()

    def run(self):
        fstd = open(self.kegg_result_path, 'w')
        ferr = open(self.kegg_error_path, 'w')
        solved = 0
        failed = 0
        for data in open(self.input_path):
            data = data.strip()
            if not data:
                continue
            try:
                if self.is_gene:
                    for output in self.work_for_gene(data):
                        fstd.write(output + '\n')
                else:
                    for output in self.work_for_kegg(data, True):
                        fstd.write(output + '\n')
                solved += 1
            except:
                ferr.write('%s\n' % data)
                failed += 1
                traceback.print_exc()
            self.logger.info(
                "Completed %d, success rate %d/%d=%.2f%%" % (
                    solved + failed, solved, solved + failed, solved * 100.0 / (solved + failed)))
            fstd.flush()
            ferr.flush()
        fstd.close()
        ferr.close()

    def work_for_gene(self, gene):
        for kegg_id in self.get_kegg_id(gene):
            for kegg_pathway in self.work_for_kegg(kegg_id, False):
                yield '%s\t%s' % (gene, kegg_pathway)

    def work_for_kegg(self, kegg_id, return_name: bool):
        names, pathways = self.get_pathway(kegg_id)
        if not return_name:
            yield '%s\t%s' % (kegg_id, '; '.join(pathways))
        else:
            for name in names:
                yield '%s\t%s\t%s' % (kegg_id, name, '; '.join(pathways))

    def get_kegg_id(self, gene):
        target_path = os.path.join(self.download_directory, 'get_kegg_id_%s.html' % gene)
        if not os.path.exists(target_path):
            url = "https://www.kegg.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&dbkey=kegg&keywords=" + gene
            for retry_time in range(3):
                try:
                    req = request.Request(url=url)
                    x = request.urlopen(req, timeout=30)
                    body = x.read()
                    for item in x.headers._headers:
                        if item[0].lower() == 'content-encoding' and item[1].lower() == 'gzip':
                            body = gzip.decompress(body)
                    body = body.decode('utf-8')
                    with open(target_path, 'w', encoding='utf8') as fw:
                        fw.write(body)
                        break
                except:
                    traceback.print_exc()
                    self.logger.info("Retry %d for %s" % (retry_time, gene))
        if not os.path.exists(target_path):
            raise ValueError("Gene not found from web: %s" % gene)
        kegg_ids = self.extract_kegg_id(target_path)
        if len(kegg_ids) == 0:
            os.remove(target_path)
            raise ValueError("Gene extract failed: %s" % gene)
        return kegg_ids

    def extract_kegg_id(self, file_path):
        with open(file_path, 'r') as f:
            body = f.read()
        parser = KeggIdHTMLParser()
        parser.feed(body)
        return parser.kegg_id_map.keys()

    def get_pathway(self, kegg_id):
        target_path = os.path.join(self.download_directory, 'get_pathway_%s.html' % kegg_id)
        if not os.path.exists(target_path):
            url = "https://www.kegg.jp/dbget-bin/www_bget?ko:" + kegg_id
            for retry_time in range(3):
                try:
                    req = request.Request(url=url)
                    x = request.urlopen(req, timeout=30)
                    body = x.read()
                    for item in x.headers._headers:
                        if item[0].lower() == 'content-encoding' and item[1].lower() == 'gzip':
                            body = gzip.decompress(body)
                    body = body.decode('utf-8')
                    with open(target_path, 'w', encoding='utf8') as fw:
                        fw.write(body)
                        break
                except:
                    traceback.print_exc()
                    self.logger.info("Retry %d for %s" % (retry_time, kegg_id))
        if not os.path.exists(target_path):
            raise ValueError("Kegg not found from web: %s" % kegg_id)
        names, pathways = self.extract_name_pathway(target_path)
        if len(pathways) == 0 or len(names) == 0:
            os.remove(target_path)
            if len(pathways) == 0:
                pathways = ["No Pathway"]
            if len(names) == 0:
                names = ["Not Found"]
        return names, pathways

    def extract_name_pathway(self, file_path):
        with open(file_path, 'r') as f:
            body = f.read()
        parser = KeggPathwayHTMLParser()
        parser.feed(body)
        return parser.names, parser.pathways
