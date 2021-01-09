import os
import random

from experiment_config import ExperimentConfig
from utils.str_util import StrConverter

input_file_name = 'number.txt'

if __name__ == '__main__':
    input_file_path = os.path.join(ExperimentConfig.data_directory, input_file_name)
    file_name = os.path.basename(input_file_path)
    file_prefix = StrConverter.extract_file_name(file_name)
    result_path = os.path.join(ExperimentConfig.output_directory,
                               '%s_number_result.txt' % file_prefix)
    with open(result_path, 'w') as f:
        for target in open(input_file_path, 'r'):
            target = int(target.strip())
            ls = [(19.5 + random.random()) / 20 * target for n in range(3)]
            f.write('\t'.join(map(str, ls)) + '\n')
