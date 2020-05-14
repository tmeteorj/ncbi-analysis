import os


class EcocycDataLoader:
    def __init__(self, ecocyc_data_file):
        self.ecocyc_data_file = ecocyc_data_file
        self.records = []
        self.inter_records = []
        self.rev_inter_records = []

    def build_database(self):
        headers = None
        inv_headers = []
        for line in open(self.ecocyc_data_file):
            items = line.rstrip('\r\n').split('\t')
            if len(items) == 0:
                continue
            if headers is None:
                for idx, header in enumerate(items):
                    headers[header] = idx
                    inv_headers.append(header)
                continue
            attr = {header: value for header, value in zip(inv_headers, items)}
            record = EcocycRecord(attr)
            direction, inter_records = record.generate_inter_record()
            self.records.append(record)
            if direction == '>':
                self.inter_records.extend(inter_records)
            else:
                self.rev_inter_records.extend(inter_records)
        self.inter_records.sort(lambda arg: arg.start)
        self.inv_inter_records.sort(lambda arg: arg.start)


class EcocycRecord:
    def __init__(self, attr):
        for header in ['gene', 'product', 'promoter_name', 'promoter_pos', 'gene_start_pos', 'map_start_pos',
                       'map_end_pos']:
            if header.endswith('pos'):
                setattr(self, header, int(attr.get(header, -1)))
            else:
                setattr(self, header, attr.get(header, ''))

    def generate_inter_record(self):
        records = []
        records.append(
            EcocycInterRecord(name=self.gene,
                              product=self.product,
                              start=self.map_start_pos,
                              end=self.map_end_pos,
                              is_gene=True)
        )
        if self.promoter_name != '':
            records.append(
                EcocycInterRecord(name=self.promoter_name,
                                  product='',
                                  start=self.promoter_pos,
                                  end=self.gene_start_pos,
                                  is_gene=False)
            )
        return records[-1].direction, records


class EcocycInterRecord:
    def __init__(self, name, product, start, end, is_gene):
        self.name = name
        self.product = product
        self.start = start
        self.end = end
        self.direction = '>' if start < end else '<'
        self.is_gene = is_gene
