import os
from dataclasses import dataclass

from utils.str_util import StrConverter


@dataclass
class LocationReorder:
    location_path: str
    index_path: str
    output_directory: str

    def __post_init__(self):
        file_name = os.path.basename(self.location_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory, '%s_reorder.txt' % file_prefix)

    def run(self):
        data = self.read_all_location()
        with open(self.result_path, 'w', encoding='utf8') as fw:
            for index in open(self.index_path, 'r', encoding='utf8'):
                index = index.strip()
                result = data.get(index, data.get('(%s)' % index, None))
                if result is None:
                    print('%s not found in location file' % index)
                for line in result:
                    fw.write(line + '\n')
                fw.write('\n')

    def read_all_location(self):
        data = {}
        buff = []
        last_index = None
        for line in open(self.location_path, 'r', encoding='utf8'):
            line = line.strip()
            if line == '': continue
            if line[0] == '(' and line[-1] == ')':
                if len(buff) > 0:
                    data[last_index] = buff[:]
                    buff.clear()
                last_index = line
            buff.append(line)
        if len(buff) > 0:
            data[last_index] = buff
        return data
