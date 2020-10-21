import os
import re
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from typing import Tuple

from src.utils.ecocyc_data_loader import EcocycDataLoader, EcocycInterRecord
from src.utils.str_util import StrConverter


@dataclass
class GeneLocationAnalysis:
    input_file_path: str
    ecocyc_file_path: str
    output_directory: str
    process_sub_data: bool = True
    filter_sub_span: Tuple[int, int] = None
    output_promoter: bool = False

    def __post_init__(self):
        self.ecocyc_data_loader = EcocycDataLoader(self.ecocyc_file_path, self.output_promoter)

        self.data_name = os.path.basename(self.input_file_path)
        file_name = os.path.basename(self.input_file_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory,
                                        '%s_location_result.txt' % file_prefix)
        self.sub_result_path = os.path.join(self.output_directory,
                                            '%s_sub_location_result.txt' % file_prefix)

    def run(self):
        self.ecocyc_data_loader.build_database()
        todo_list = []
        buff = []
        for line in open(self.input_file_path, 'r', encoding='utf8'):
            line = line.strip()
            if line == '':
                continue
            elif line[0] == '(' and line[-1] == ')':
                continue
            elif line.startswith('>NC'):
                if len(buff) > 0:
                    todo_list.append(self.parse_similarity_data(buff))
                    buff.clear()
            buff.append(line)
        if len(buff) > 0:
            todo_list.append(self.parse_similarity_data(buff))
        with open(self.result_path, 'w', encoding='utf8') as fw:
            for idx, data in enumerate(todo_list):
                self.process_one_data(data)
                fw.write('(%d)\n' % (idx + 1))
                fw.write(data['header'] + '\n')
                fw.write(data['match_info'] + '\n')
                fw.write(data['direction'] + '\n')
                for line in data['additional']:
                    fw.write(line + '\n')
                for location_info in data['location_result']:
                    fw.write(location_info + '\n')
                fw.write('\n')
        if self.process_sub_data:
            with open(self.sub_result_path, 'w', encoding='utf8') as fw:
                for idx, data in enumerate(todo_list):
                    for sub_idx, sub_data in enumerate(self.extract_sub_data(data)):
                        self.process_one_data(sub_data)
                        if self.pass_filter_sub_location(sub_data):
                            fw.write('(%d-%d)\n' % (idx + 1, sub_idx + 1))
                            fw.write(sub_data['header'] + '\n')
                            fw.write('Original Position\t' + sub_data['left'] + '\t' + sub_data['right'] + '\n')
                            fw.write(sub_data['match_info'] + '\n')
                            fw.write(sub_data['direction'] + '\n')
                            for line in sub_data['additional']:
                                fw.write(line + '\n')
                            for location_info in sub_data['location_result']:
                                fw.write(location_info + '\n')
                            fw.write('\n')

    def pass_filter_sub_location(self, sub_data):
        def extract_gene_start_end(gene_sequence):
            items = re.split(r'\t|\n', gene_sequence)
            for item in items:
                if re.match('\d+-\d+', item):
                    return tuple(map(int, item.split('-')))
            return None, None

        if not self.filter_sub_span:
            return True
        if sub_data['location_result'][0].startswith('5\''):
            gene_start, gene_end = extract_gene_start_end(sub_data['location_result'][0])
            if gene_start < gene_end:
                span_start, span_end = gene_start - self.filter_sub_span[0], gene_start + self.filter_sub_span[1]
                return span_start <= sub_data['start'] <= sub_data['end'] <= span_end
            else:
                span_start, span_end = gene_start + self.filter_sub_span[0], gene_start - self.filter_sub_span[1]
                return span_start >= sub_data['start'] >= sub_data['end'] >= span_end
        return True

    def extract_sub_data(self, data):
        match_info = data['match_info'].split('\n')
        best_cnt = None
        for kv in match_info:
            if kv.find(':') >= 0:
                k, v = kv.split(':')
                if k.find('consistency_match_format') >= 0:
                    match_format = v.strip()
            elif kv.find('consistency_similarity') >= 0:
                k, v = kv.split('\t')
                best_cnt = int(float(v.strip()))
        cur_cnt = 0
        start = None
        step = 1 if data['start'] < data['end'] else -1
        for end, m in enumerate(match_format):
            if m == '*':
                if cur_cnt == 0:
                    start = end
                cur_cnt += 1
            elif m == '.':
                cur_cnt = 0
            if cur_cnt == best_cnt:
                sub_start = data['start'] + step * start
                sub_end = data['start'] + step * end
                sub_data = deepcopy(data)
                sub_data['left'] = str(start + 1)
                sub_data['right'] = str(end + 1)
                sub_data['start'] = sub_start
                sub_data['end'] = sub_end
                sub_data['header'] = '%s/%s-%s' % (data['header'].split('/')[0], sub_start, sub_end)
                output = []
                for match_info_data in match_info:
                    if match_info_data.find(':') >= 0:
                        k, v = match_info_data.split(':')
                        if k.find('consistency') >= 0:
                            output.append(k + ': ' + v.strip()[start:end + 1])
                        else:
                            output.append(k + ': ' + v.strip())
                    else:
                        k, v = match_info_data.split('\t')
                        output.append(k + '\t' + v)
                sub_data['match_info'] = '\n'.join(output)
                sub_data['location_result'].clear()
                yield sub_data

    def process_one_data(self, data):
        start, end = data['start'], data['end']
        idx = self.ecocyc_data_loader.find_first_le(start)
        data['location_result'].extend(
            self.get_location_information(self.ecocyc_data_loader.inter_records, idx, start, end))

    def get_location_information(self, inter_records, idx, start, end):
        left = min(start, end)
        right = max(start, end)
        find_left = max(idx - 2, 0)
        find_right = min(idx + 3, len(inter_records))
        result = []
        left_neareast_record = None
        right_neareast_record = None
        for index in range(find_left, find_right):
            intersect_status = ''
            record = inter_records[index]
            status = interval_check(record.left, record.right, left, right)
            if status in [IntervalPositionStatus.IntersectLeft, IntervalPositionStatus.CoverLeft]:
                if record.direction == '>':
                    intersect_status = '5\''
                else:
                    intersect_status = '3\''
            elif status in [IntervalPositionStatus.IntersectRight, IntervalPositionStatus.CoverRight]:
                if record.direction == '>':
                    intersect_status = '3\''
                else:
                    intersect_status = '5\''
            elif status in [IntervalPositionStatus.Inner]:
                intersect_status = 'cds'
            elif status in [IntervalPositionStatus.Cover]:
                intersect_status = 'cover'
            elif status in [IntervalPositionStatus.TotallyLeft]:
                intersect_status = 'inter-genic'
                if right_neareast_record is None or right_neareast_record.left > record.left:
                    right_neareast_record = record
            elif status in [IntervalPositionStatus.TotallyRight]:
                intersect_status = 'inter-genic'
                if left_neareast_record is None or left_neareast_record.right < record.right:
                    left_neareast_record = record
            else:
                raise ValueError('status not correct')

            if intersect_status != 'inter-genic':
                result.append(self.render_location_result(intersect_status, record, left, right))
        left_name = 'None' if not left_neareast_record else left_neareast_record.name
        right_name = 'None' if not right_neareast_record else right_neareast_record.name
        if len(result) == 0:
            result.append('inter-genic of %s, %s' % (left_name, right_name))
        assert len(result) > 0
        return result

    @staticmethod
    def render_location_result_inter_genic(left_record, right_record):
        result = ['inter-genic of %s, %s' % (left_record.name, right_record.name)]
        for record in [left_record, right_record]:
            result.append(
                'direction of %s[%d-%d]: %s\n' % (record.name, record.start, record.end, record.direction * 10))

    @staticmethod
    def render_location_result(intersect_status, record: EcocycInterRecord, left, right):
        result = ['%s of %s' % (intersect_status, record.name)]
        intersect_directions = ''
        original_directions = ''
        record_len = record.right - record.left + 1
        part_len = record_len // 10
        max_draw = 10
        if part_len == 0:
            part_len = 1
            max_draw = record_len
        direction = record.direction
        for idx in range(max_draw):
            r_left = record.left + idx * part_len
            r_right = r_left + part_len - 1 if idx < max_draw - 1 else record.right
            coverage = count_coverage([left, right], [r_left, r_right]) * 100.0 / part_len
            if coverage > 50.0:
                intersect_directions += '*'
            else:
                intersect_directions += direction
            original_directions += direction
        result.append('original direction  : %s' % original_directions)
        result.append('intersect direction : %s' % intersect_directions)
        if record.is_gene:
            result.append('%s-%s\tgene=%s\tproduct=%s' % (record.start, record.end, record.name, record.product))
        else:
            result.append('%s-%s\tpromoter=%s' % (record.start, record.end, record.name))
        return '\n'.join(result)

    @staticmethod
    def parse_similarity_data(buff):
        data = {
            'additional': buff[1:],
        }
        primary_info, match_info = buff[0].split('\t')
        file_info, inter_info = primary_info.split('/')
        start, end = inter_info.split('-')
        data['start'] = int(start)
        data['end'] = int(end)
        data['header'] = '%s/%s-%s' % (file_info, start, end)
        match_info = re.sub(r'direction=(\+|\-),', '', match_info).split(',')
        match_info = [kv.split('=') for kv in match_info]
        output = []
        for k, v in match_info:
            if k.find(':') >= 0:
                output.append(k + v)
            else:
                output.append(k + '\t' + v)
        data['match_info'] = '\n'.join(output)
        data['direction'] = ('>' if start < end else '<') * 10
        data['location_result'] = []
        return data


def count_coverage(seg_a, seg_b):
    if seg_a[0] > seg_b[0]:
        seg_a, seg_b = seg_b, seg_a
    if seg_b[1] <= seg_a[1]:
        return seg_b[1] - seg_b[0] + 1
    elif seg_b[0] <= seg_a[1]:
        return seg_a[1] - seg_b[0] + 1
    else:
        return 0


def interval_check(record_left, record_right, left, right):
    if right < record_left:
        return IntervalPositionStatus.TotallyLeft

    elif left < record_left <= right < record_right:
        return IntervalPositionStatus.IntersectLeft

    elif left < record_left <= record_right <= right:
        return IntervalPositionStatus.CoverLeft

    elif record_left <= left <= right <= record_right:
        return IntervalPositionStatus.Inner

    elif left <= record_left <= record_right < right:
        return IntervalPositionStatus.CoverRight

    elif record_left < left <= record_right < right:
        return IntervalPositionStatus.IntersectRight

    elif record_right < left:
        return IntervalPositionStatus.TotallyRight

    elif left < record_left <= record_right < right:
        return IntervalPositionStatus.Cover

    else:
        raise ValueError("[%d,%d] <-> [%d,%d]" % (record_left, record_right, left, right))


def format_data_to_tsv(input_path, output_path, ecocyc_data_loader):
    headers = ['index', 'weighted_similarity', 'text_distance_similarity',
               'direct_match_similarity', 'consistency_similarity', 'location',
               'gene_name', 'type', 'exonic_gene_sizes', 'product',
               'site']
    max_header_len = len(headers)
    datas = []
    buff = []
    for line in open(input_path, 'r', encoding='utf8'):
        line = line.strip()
        if line == '':
            if len(buff) > 0:
                for data in extract_consistency_record(buff, ecocyc_data_loader):
                    if data is not None:
                        output = []
                        for header in headers:
                            output.append(data.get(header, ''))
                        for start, end in data.get('location interval', []):
                            output.extend([start, end])
                        max_header_len = max(max_header_len, len(output))
                        datas.append(output)
                buff.clear()
            continue
        buff.append(line)
    if len(buff) > 0:
        for data in extract_consistency_record(buff, ecocyc_data_loader):
            if data is not None:
                output = []
                for header in headers:
                    output.append(data.get(header, ''))
                max_header_len = max(max_header_len, len(output))
                datas.append(output)
    with open(output_path, 'w', encoding='utf8') as fw:
        idx = 1
        while len(headers) < max_header_len:
            headers.extend(['location start %d' % idx, 'location end %d' % idx])
            idx += 1
        fw.write('\t'.join(headers) + '\n')
        for data in datas:
            fw.write('\t'.join(data) + '\n')


def extract_consistency_record(buff, ecocyc_data_loader: EcocycDataLoader):
    def update_data(odata, location_type, genes, direction_matched, direction):
        data = {k: v for k, v in odata.items()}
        if location_type == 'inter-genic':
            data['location'] = 'inter genic'
            data['gene_name'] = genes
        else:
            data['location'] = 'antisense' if direction_matched == direction else 'sense'
            if location_type == '5\'' or location_type == '3\'':
                data['location'] += ' ' + location_type + 'utr'
            else:
                data['location'] += ' ' + location_type
            data['gene_name'] = genes
            record = ecocyc_data_loader.get_target_gene(genes.strip())
            if record is None:
                print(genes + ' not found, might be a promoter')
            else:
                data['type'] = record.type
                data['exonic_gene_sizes'] = record.exonic_gene_sizes
                data['product'] = record.product
        return data

    data = {}
    location_type = ''
    direction = None
    direction_matched = None
    genes = ''
    for line in buff:
        items = line.split('\t')
        if items[0] in ['weighted_similarity',
                        'text_distance_similarity',
                        'direct_match_similarity',
                        'consistency_similarity']:
            data[items[0]] = line.split('\t')[1].strip('%')
        elif line.startswith('>>>'):
            direction = '>'
        elif line.startswith('<<<'):
            direction = '<'
        elif line.find(' of ') >= 0:
            items = line.split(' of ')
            if len(items) != 2 or items[0] not in ['5\'', '3\'', 'cds', 'cover', 'inter-genic']:
                continue
            if location_type != '':
                yield update_data(data, location_type, genes, direction_matched, direction)
            location_type = items[0]
            genes = items[1]
        elif line.startswith('original direction'):
            direction_matched = line[-1]
        elif line.startswith('>NC'):
            data['site'] = line.strip().split('/')[-1]
        elif line.startswith('(') and line.strip().endswith(')'):
            data['index'] = line.strip()[1:-1]
        elif line.startswith('match_format'):
            _, sequence = re.split(r'\s+', line.strip())
            sequence = sequence[1:]
            score = 0
            interval = []
            for idx, ch in enumerate(sequence):
                if ch == '.':
                    score = 0
                else:
                    score += 1
                if score == int(data.get('consistency', 10000)):
                    interval.append([str(idx + 2 - score), str(idx + 1)])
            data['location interval'] = interval
    yield update_data(data, location_type, genes, direction_matched, direction)


class IntervalPositionStatus(Enum):
    # .. --
    TotallyLeft = 0
    # ..*--
    IntersectLeft = 1
    # ..***
    CoverLeft = 2
    # -***-
    Inner = 3
    # ***..
    CoverRight = 4
    # -**..
    IntersectRight = 5
    # -- ..
    TotallyRight = 6
    # .***.
    Cover = 7
