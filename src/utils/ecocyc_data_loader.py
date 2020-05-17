class EcocycDataLoader:
    def __init__(self, ecocyc_data_file):
        self.ecocyc_data_file = ecocyc_data_file
        self.gene_name_index = {}
        self.records = []
        self.inter_records = []
        self.headers = None
        self.inv_headers = []

    def build_database(self):
        for line in open(self.ecocyc_data_file, 'r', encoding='utf8', errors='ignore'):
            items = line.rstrip('\r\n').split('\t')
            if len(items) == 0:
                continue
            if self.headers is None:
                self.headers = {}
                for idx, header in enumerate(items):
                    self.headers[header] = idx
                    self.inv_headers.append(header)
                continue
            attr = {header: value for header, value in zip(self.inv_headers, items)}
            record = EcocycRecord(attr)
            direction, inter_records = record.generate_inter_record()
            self.records.append(record)
            self.inter_records.extend(inter_records)
            self.gene_name_index[record.gene] = len(self.records) - 1
        self.inter_records.sort(key=lambda arg: arg.start)

    def find_first_le(self, pos):
        return binary_search_first_le(self.inter_records,
                                      0,
                                      len(self.inter_records) - 1,
                                      pos)

    def get_target_gene(self, name):
        try:
            return self.records[self.gene_name_index[name]]
        except:
            return None


def binary_search_first_le(arr, left, right, value):
    while left < right:
        mid = (left + right) // 2
        if arr[mid].start >= value:
            right = mid
        elif arr[mid].start < value:
            left = mid + 1
    return left


class EcocycRecord:
    def __init__(self, attr):
        for header in ['gene', 'product_type', 'product', 'promoter_name', 'promoter_pos', 'gene_start_pos',
                       'map_start_pos',
                       'map_end_pos', 'exonic_gene_sizes', 'type']:
            if header.endswith('pos'):
                value = attr.get(header, -1)
                if value == '':
                    value = -1
                setattr(self, header, int(value))
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

        self.left = min(self.start, self.end)
        self.right = max(self.start, self.end)
