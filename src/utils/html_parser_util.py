import logging
import re
from html.parser import HTMLParser

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

ecocyc_id_script_prefix = 'window.location.replace(\'/gene?'
inf = 1000000


class BaseHTMLParser(HTMLParser):
    @staticmethod
    def extract_map_position(data):
        start = data.index('[')
        end = data.index(']')
        data = data[start + 1:end]
        if data.find('<-') > 0:
            end, start = data.split('<-')
        else:
            start, end = data.split('->')
        start, end = int(start.replace(',', '')), int(end.replace(',', ''))
        return start, end

    @staticmethod
    def extract_gene_name(data):
        return re.sub(r'<\w+>', '', data)

    @staticmethod
    def extract_id_from_script(data: str):
        if data.find('gene:\'') < 0:
            return None
        start = data.index('gene:\'') + 6
        end = data.index('\'', start)
        return data[start:end]

    @staticmethod
    def extract_id_from_data(data):
        items = re.split(r'\'|\?|&|\"', data)
        for kv in items:
            if kv.find('=') > 0:
                k, v = kv.split('=', 1)
                if k == 'id':
                    return v
        return None


class EcocycHTMLParser(BaseHTMLParser):
    def __init__(self, do_extract_id=False, gene_name=None, do_extract_summary=False):
        super(EcocycHTMLParser, self).__init__()
        self.last_td_data = None
        self.last_ecocyc_id = None
        self.last_a_data = None

        self.do_extract_id = do_extract_id
        self.do_extract_summary = do_extract_summary
        self.gene_name = gene_name
        self.depth = 0
        self.fill_depth = -inf

        self.extract_attr = {'location': None, 'reaction': None, 'gene': None, 'enzyme': None, 'rna': None,
                             'protein': None, 'polypeptide': None, 'function when intact': None, 'transporter': None,
                             'map position': None, 'summary': None}
        self.ecocyc_id = None
        self.summary_extract_step = 'not_start'

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.last_a_data = ''
        if tag == 'td':
            self.depth += 1
        tag = tag.strip()
        logger.debug("Start tag: %s" % tag)
        if tag == 'a' and self.do_extract_id:
            href = None
            for attr in attrs:
                if attr[0] == 'href':
                    href = attr[1]
            if href is not None and href.startswith('/gene?orgid=ECOLI&id='):
                self.last_ecocyc_id = self.extract_id_from_data(href)
        if tag == 'p' and self.do_extract_summary and self.summary_extract_step == 'start':
            for attr in attrs:
                if attr[0] == 'class' and attr[1] == 'ecoparagraph':
                    self.summary_extract_step = 'end'

    def handle_endtag(self, tag):
        if tag == 'td':
            self.depth -= 1
            if self.depth < self.fill_depth:
                if self.extract_attr[self.last_td_data] != '':
                    if self.last_td_data == 'map position':
                        map_data = self.extract_attr[self.last_td_data]
                        self.extract_attr[self.last_td_data] = self.extract_map_position(map_data)
                    self.fill_depth = -inf
                    self.last_td_data = None
        elif tag == 'a':
            if self.last_td_data == 'reaction' and self.extract_attr['reaction']:
                self.extract_attr['reaction'] += '__#####__'
            if self.do_extract_id and self.last_ecocyc_id is not None:
                gene_name = self.extract_gene_name(self.last_a_data)
                if gene_name == self.gene_name:
                    self.ecocyc_id = self.last_ecocyc_id
                else:
                    self.last_ecocyc_id = None
            self.last_a_data = None
        logger.debug("End tag  :%s" % tag)

    def handle_data(self, data):
        data = data.strip()
        if self.last_a_data is not None:
            self.last_a_data += data
        if data == 'Locations' or data == 'Reactions':
            data = data[:-1]
        if self.do_extract_summary:
            if self.summary_extract_step == 'not_start' and data == 'Summary' and self.lasttag == 'h3':
                self.extract_attr['summary'] = ''
                self.summary_extract_step = 'start'
            elif self.summary_extract_step == 'start':
                self.extract_attr['summary'] += data
        elif self.do_extract_id:
            if self.lasttag == 'script' and data.startswith(ecocyc_id_script_prefix):
                data = data[len(ecocyc_id_script_prefix):]
                self.ecocyc_id = self.extract_id_from_data(data)
        elif data != '':
            if self.last_td_data in self.extract_attr:
                self.extract_attr[self.last_td_data] += data
            if self.lasttag == 'td':
                if self.fill_depth == -inf:
                    self.last_td_data = data.lower()
                    if self.last_td_data in self.extract_attr:
                        self.fill_depth = self.depth
                        self.extract_attr[self.last_td_data] = ''
            if data.find('typeObjectPage') > 0:
                self.ecocyc_id = self.extract_id_from_script(data)
            logger.debug("Data     :%s" % data)


class UrlHTMLParser(BaseHTMLParser):
    def __init__(self):
        super(UrlHTMLParser, self).__init__()
        self.ecocycs = []

    def handle_starttag(self, tag, attrs):
        tag = tag.strip()
        logger.debug("Start tag: %s" % tag)
        if tag == 'a':
            href = None
            for attr in attrs:
                if attr[0] == 'href':
                    href = attr[1]
            if href is not None:
                href = href.replace('&amp;', '&')
                self.ecocycs.append([href + '#tab=TU', self.extract_name_from_data(href), ''])

    def handle_endtag(self, tag):
        logger.debug("End tag  :%s" % tag)

    def handle_data(self, data):
        data = data.strip()
        if data != '':
            self.ecocycs[-1][-1] += data

    @staticmethod
    def extract_name_from_data(data):
        items = re.split(r'\'|\?|&|\"', data)
        for kv in items:
            if kv.find('=') > 0:
                k, v = kv.split('=', 1)
                if k in ['id', 'object']:
                    return v
        return None


class GoHTMLParser(BaseHTMLParser):
    def __init__(self):
        super(GoHTMLParser, self).__init__()

        self.tb_depth = 0
        self.tag_stack = []
        self.tr_depth = []
        self.td_depth = []
        self.go_table = []

    def handle_starttag(self, tag, attrs):
        self.tag_stack.append(tag)
        tag = tag.strip()
        if tag == 'table':
            self.tb_depth += 1
            self.tr_depth.append(0)
            self.td_depth.append(0)
            if self.tb_depth == 1 and len(list(filter(lambda arg: arg[0] == 'class', attrs))) == 0:
                self.tb_depth = 1000
        elif tag == 'td':
            self.td_depth[-1] += 1
        elif tag == 'tr':
            self.tr_depth[-1] += 1
            self.td_depth[-1] = 0
        logger.debug("Start tag  :%s" % tag)

    def handle_endtag(self, tag):
        self.tag_stack.pop(-1)
        if tag == 'table':
            self.tb_depth -= 1
            self.tr_depth.pop(-1)
            self.td_depth.pop(-1)
            if self.tb_depth == 0:
                self.tb_depth = 1000
        logger.debug("End tag  :%s" % tag)

    def handle_data(self, data):
        if self.tb_depth == 1 and self.td_depth[-1] == 1 and self.tag_stack[-1] == 'td':
            data = re.sub(r'^\s+', '', data)
            data = re.sub(r'(\s|:)+$', '', data)
            self.go_table.append([data, ''])
        elif self.tb_depth == 2 and self.td_depth[-1] == 2 and self.tag_stack[-1] == 'a':
            self.go_table[-1][-1] = (self.go_table[-1][-1] + ',' + data.strip()).lstrip(',')
        logger.debug("Data     :%s" % data)
