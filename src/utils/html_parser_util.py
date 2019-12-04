import logging
import re
from html.parser import HTMLParser

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

ecocyc_id_script_prefix = 'window.location.replace(\'/gene?'
inf = 1000000


class EcocycHTMLParser(HTMLParser):
    def __init__(self, do_extract_id=False, gene_name=None):
        super(EcocycHTMLParser, self).__init__()
        self.last_td_data = None
        self.last_ecocyc_id = None

        self.do_extract_id = do_extract_id
        self.gene_name = gene_name
        self.depth = 0
        self.fill_depth = -inf

        self.extract_attr = {'Location': None, 'Reaction': None}
        self.ecocyc_id = None

    def handle_starttag(self, tag, attrs):
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

    def handle_endtag(self, tag):
        if tag == 'td':
            self.depth -= 1
            if self.depth < self.fill_depth:
                if self.extract_attr[self.last_td_data] is None:
                    self.extract_attr[self.last_td_data] = ''
                else:
                    self.fill_depth = -inf
                    self.last_td_data = None

        logger.debug("End tag  :%s" % tag)

    def handle_data(self, data):
        data = data.strip()
        if data == 'Locations' or data == 'Reactions':
            data = data[:-1]
        if self.do_extract_id:
            if self.lasttag == 'script' and data.startswith(ecocyc_id_script_prefix):
                data = data[len(ecocyc_id_script_prefix):]
                self.ecocyc_id = self.extract_id_from_data(data)
            elif self.lasttag in ['a', 'b'] and self.last_ecocyc_id is not None:
                gene_name = self.extract_gene_name(data)
                if gene_name == self.gene_name:
                    self.ecocyc_id = self.last_ecocyc_id
                else:
                    self.last_ecocyc_id = None
        elif data != '':
            if self.last_td_data in self.extract_attr:
                self.extract_attr[self.last_td_data] += data
            if self.lasttag == 'td':
                if self.fill_depth == -inf:
                    self.last_td_data = data
                    if data in self.extract_attr:
                        self.fill_depth = self.depth
            if data.find('typeObjectPage') > 0:
                self.ecocyc_id = self.extract_id_from_script(data)
            logger.debug("Data     :%s" % data)

    @staticmethod
    def extract_gene_name(data):
        return re.sub(r'<\w+>', '', data)

    @staticmethod
    def extract_id_from_script(data: str):
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


class UrlHTMLParser(HTMLParser):
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
