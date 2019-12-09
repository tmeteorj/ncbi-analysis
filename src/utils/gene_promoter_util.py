import re


class GeneTUInfo:
    def __init__(self, items=None):
        self.idx = None
        self.link = None
        self.left = None
        self.top = None
        self.right = None
        self.bottom = None
        self.attribute = None
        if items is not None:
            self.parse_from_list(items)

    def parse_from_list(self, items):
        assert len(items) == 7, 'items size in table unit not correct'
        self.idx = int(items[0])
        self.link = items[1]
        self.left, self.top, self.right, self.bottom = map(int, items[2:6])
        self.attribute = self.parse_attributes(items[6])

    def is_gene(self, gene_name=None):
        if gene_name is not None:
            gene_attr = self.attribute.get('Gene', '').split()
            return gene_name in gene_attr
        else:
            return 'Gene' in self.attribute

    def is_promoter(self, check_start_site=False):
        if check_start_site:
            return 'Promoter' in self.attribute and 'Tr.Start site' in self.attribute
        else:
            return 'Promoter' in self.attribute

    def get_promoter_name(self):
        return self.attribute.get('Promoter', None)

    def get_promoter_start_site(self, int_pos=False):
        position = self.attribute.get('Tr.Start site', None)
        if int_pos:
            position = int(position.replace(',', ''))
        return position

    def get_gene_start_position(self):
        a, op, b = re.split(r'\s+', self.attribute['Location'])
        if op == '<-':
            return int(b.replace(',', ''))
        else:
            return int(a.replace(',', ''))

    def get_direction_of_gene(self):
        location = self.attribute['Location']
        if location.find('<-') > 0:
            return 'Left'
        elif location.find('->') > 0:
            return 'Right'
        else:
            raise RuntimeError('Get direction of gene failed')

    def same(self, x):
        return x.idx == self.idx

    def __str__(self):
        result = '\t'.join(map(str, [self.left, self.top, self.right, self.bottom]))
        if self.is_promoter():
            result += '\t' + self.attribute['Promoter']
        else:
            result += 'Other'
        return result

    @staticmethod
    def parse_attributes(attr_str):
        attr = re.sub(r'<b>|</b>', '', attr_str)
        result = {}
        for line in re.split('<BR>|<br>', attr):
            try:
                if line.find(':') < 0:
                    continue
                k, v = map(lambda arg: arg.strip(), line.split(':', 1))
                result[k] = v
            except:
                print('Parse gene table unit error for ' + line)
        if len(result) == 0:
            result['Body'] = attr_str
        return result


def get_all_promoters(tu_list, check_start_site=False):
    return list(filter(lambda arg: arg.is_promoter(check_start_site), tu_list))


def get_all_genes(tu_list, direction=None):
    return list(
        filter(lambda arg: arg.is_gene() and (direction is None or direction == arg.get_direction_of_gene()), tu_list))


def get_valid_range(target_gene: GeneTUInfo, gene_list: list):
    direction = target_gene.get_direction_of_gene()
    if direction == 'Left':
        bigger_than = target_gene.right
        less_than = None
        for gene_tu in gene_list:
            if direction == gene_tu.get_direction_of_gene():
                if gene_tu.left > bigger_than:
                    if less_than is None or less_than > gene_tu.left:
                        less_than = gene_tu.left
        return 'Left', bigger_than, less_than
    else:
        less_than = target_gene.left
        bigger_than = None
        for gene_tu in gene_list:
            if direction == gene_tu.get_direction_of_gene():
                if gene_tu.right < less_than:
                    if bigger_than is None or bigger_than < gene_tu.right:
                        bigger_than = gene_tu.right
        return 'Right', bigger_than, less_than


def filter_same_direction(gene_tu: GeneTUInfo, tu_list: list):
    top = gene_tu.top
    result = []
    for tu in tu_list:
        if tu.top <= top:
            result.append(tu)
    return result


def get_better_one(promoter_a, promoter_b, direction):
    position_a = promoter_a.get_promoter_start_site(True)
    position_b = promoter_b.get_promoter_start_site(True)
    comp = position_a < position_b
    if direction == 'Left':
        return promoter_b if comp else promoter_a
    elif direction == 'Right':
        return promoter_a if comp else promoter_b
    else:
        raise ValueError('Direction error')


def get_target_promoter(target_gene: GeneTUInfo, tu_list: list):
    direction = target_gene.get_direction_of_gene()
    genes = get_all_genes(tu_list, direction)
    promoters = get_all_promoters(tu_list, True)
    promoters = filter_same_direction(target_gene, promoters)
    ls = genes + promoters
    ls.sort(key=lambda arg: (arg.left if direction == 'Right' else arg.right) * 10 + int(arg.is_gene()))
    tot = len(ls)
    add = 1 if direction == 'Right' else -1
    idx = 0 if direction == 'Right' else (tot - 1)
    last_promoter = None
    gene_appears = False
    while 0 <= idx < tot:
        item = ls[idx]
        if item.is_gene():
            gene_appears = True
            if item.same(target_gene):
                return last_promoter
        elif item.is_promoter():
            if gene_appears or last_promoter is None:
                last_promoter = item
                gene_appears = False
        idx += add
    return None

#
# def get_target_promoter(target_gene: GeneTUInfo, tu_list: list):
#     genes = get_all_genes(tu_list)
#     direction, bigger_than, less_than = get_valid_range(target_gene, genes)
#     if bigger_than is None: bigger_than = -1000000000000
#     if less_than is None: less_than = 10000000000000
#     promoters = get_all_promoters(tu_list, direction, bigger_than, less_than, check_start_site=True)
#     promoters = filter_same_direction(target_gene, promoters)
#     best = None
#     for promoter in promoters:
#         if best is None:
#             best = promoter
#         else:
#             best = get_better_one(best, promoter, direction)
#     return best
