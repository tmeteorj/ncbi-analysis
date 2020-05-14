import os
from enum import Enum

from src.utils.ecocyc_data_loader import EcocycDataLoader
from src.utils.str_util import StrConverter


class GeneLocationAnalysis:

    def __init__(self, input_file_path, ecocyc_file_path, outut_directory):
        self.input_file_path = input_file_path
        self.output_directory = outut_directory
        self.ecocyc_data_loader = EcocycDataLoader(ecocyc_file_path)

        self.data_name = os.path.basename(input_file_path)
        file_name = os.path.basename(input_file_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory,
                                        '%s_location_result.txt' % file_prefix)

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

    def process_one_data(self, data):
        start, end = data['start'], data['end']
        idx, rev_idx = self.ecocyc_data_loader.find_first_le(start)

    def get_location_information(self, inter_records, idx, start, end):
        if inter_records[idx].direction == '>' and start < end:
            same_direction = True
        elif inter_records[idx].direction == '<' and start > end:
            same_direction = True
        else:
            same_direction = False

        left = min(start, end)
        right = max(start, end)
        find_left = max(idx - 2, 0)
        find_right = min(idx + 3, len(inter_records))
        for index in range(find_left, find_right):
            record = inter_records[index]
            status = interval_check(record.left, record.right, left, right)

    @staticmethod
    def parse_similarity_data(buff):
        data = {
            'additional': buff[1:],
        }
        primary_info, match_info = buff[0].split('\t')
        file_info, inter_info = primary_info.split('/')
        end, start = inter_info.split('-')
        data['start'] = int(start)
        data['end'] = int(end)
        data['header'] = '%s/%s-%s' % (file_info, start, end)
        data['match_info'] = match_info
        return data


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
