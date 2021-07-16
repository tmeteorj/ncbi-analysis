import os
import random
import re
import time
import traceback
from collections import Counter
from dataclasses import dataclass

from experiment_config import ExperimentConfig
from utils.data_download_util import DataDownloadTool
from utils.factories.logger_factory import LoggerFactory
from utils.ncbi_database import NCBIDatabase
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter


@dataclass
class NeighborAnalysis:
    input_path: str
    download_directory: str
    output_directory: str
    keep_prefix_num: int = 1

    def __post_init__(self):
        self.logger = LoggerFactory(3)

        file_name = os.path.basename(self.input_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.neighbor_result_path = os.path.join(self.output_directory, '%s_neighbor_result.txt' % file_prefix)
        self.next_gene_result_path = os.path.join(self.output_directory, '%s_next_neighbor_result.txt' % file_prefix)
        self.source_count_path = os.path.join(self.output_directory, '%s_source_count_result.txt' % file_prefix)
        self.gene_count_path = os.path.join(self.output_directory, '%s_gene_count_result.txt' % file_prefix)

        error_directory = os.path.join(self.output_directory, 'error')
        if not os.path.exists(error_directory):
            os.makedirs(error_directory)
        self.error_result_path_prefix = os.path.join(error_directory, '%s_error_result' % file_prefix)

    def download_and_analysis(self, key, inter, file_path):
        self.logger.info('\nstart working for ' + key)
        try:
            if not os.path.exists(file_path):
                if not DataDownloadTool.download_data(key, file_path):
                    return False, None
            flag, data = self.analysis_download_file(file_path, inter)
            if flag:
                return True, data
        except:
            traceback.print_exc()
            return False, None
        return False, None

    def find_neighbor_batch(self, datas, iteration_time):
        fw = open(self.neighbor_result_path, 'a')
        solve_cnt, success_cnt, total_cnt = 0, 0, len(datas)
        logger = LoggerFactory(1)
        logger.info_with_expire_time(
            '[Iteration %d]completed %d/%d=%.2f%%' % (
                iteration_time, solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt),
            solve_cnt, total_cnt)
        fe = open(self.error_result_path_prefix + ".iter-%d.txt" % iteration_time, 'w')
        fail_datas = []
        for key, inter, additional in datas:
            solve_cnt += 1
            file_path = os.path.join(self.download_directory, key + '.txt')
            flag, data = self.download_and_analysis(key, inter, file_path)
            if flag:
                success_cnt += 1
                direction = '+' if (inter[0] < inter[1]) else '-'
                fw.write('>%s/%s-%s(%s)\n' % (key, inter[0], inter[1], direction))
                if additional != '':
                    for kv in additional.split(','):
                        k, v = kv.split('=')
                        fw.write('%s\t%s\n' % (k, v))
                fw.write('SOURCE\t%s\n' % (data.get('source', 'UNKNOWN')))
                for elem in data['data']:
                    fw.write('%s\n' % elem)
                fw.write('sequence\t%s\n' % (data.get('sequence', '')))
                fw.write('\n')
                fw.flush()
            else:
                if os.path.exists(file_path):
                    os.remove(file_path)
                fe.write('>%s/%s-%s\n' % (key, inter[0], inter[1]))
                fe.flush()
                fail_datas.append([key, inter])
            self.logger.info_with_expire_time('[Iteration %d]completed %d/%d=%.2f%%, success %d/%d=%.2f%%' % (
                iteration_time, solve_cnt, total_cnt, solve_cnt * 100.0 / total_cnt, success_cnt, solve_cnt,
                success_cnt * 100.0 / solve_cnt),
                                              solve_cnt, total_cnt)
            time.sleep(random.random())
        self.logger.info('[Iteration %d]done .' % iteration_time)
        fw.close()
        return fail_datas

    def extract_data(self, buff):
        name = buff[0][1:].strip()
        name, inter = name.split('/')
        direction = inter[-2]
        left, right = map(int, inter[:-3].split('-'))
        _, source = buff[1].strip().split('\t')
        source = self.get_prefix(source)
        target = None
        for line in buff[2:]:
            try:
                gene = self.read_gene(line)
                if self.check_gene(left, right, direction, gene, target):
                    target = gene
            except:
                continue
        return {
            'name': name,
            'direction': direction,
            'left': left,
            'right': right,
            'source': source,
            'gene': target
        }

    def get_prefix(self, source):
        if self.keep_prefix_num > 0:
            return ' '.join(re.split('\s+', source)[:self.keep_prefix_num])
        return source

    def source_gene_distribution_analysis(self):
        self.logger.info('Start source_gene_distribution_analysis')
        datas = []
        buff = []
        for line in open(self.neighbor_result_path, 'r'):
            if len(line.strip()) == 0:
                if len(buff) > 0:
                    datas.append(self.extract_data(buff))
                    buff = []
            else:
                buff.append(line.strip())
        if len(buff) > 0:
            datas.append(self.extract_data(buff))
        source_counter = Counter()
        gene_counter = Counter()
        with open(self.next_gene_result_path, 'w') as fw:
            for data in datas:
                if data['gene'] is None:
                    continue
                fw.write('>%s/%s-%s(%s)\n' % (data['name'], data['left'], data['right'], data['direction']))
                fw.write('SOURCE\t%s\n' % (data['source']))
                fw.write('%s-%s\t%s\n\n' % (data['gene']['left'], data['gene']['right'], data['gene']['gene']))
                source_counter[data['source']] += 1
                gene_counter[data['gene']['gene']] += 1
        total = len(datas)
        for file_path, counter in [(self.source_count_path, source_counter), (self.gene_count_path, gene_counter)]:
            with open(file_path, 'w') as fw:
                for k, v in counter.most_common():
                    fw.write('%s\t%d\t%.4f%%\n' % (k, v, v * 100.0 / total))

        self.logger.info('End source_gene_distribution_analysis')

    def run(self):
        open(self.neighbor_result_path, 'w').close()
        with open(self.input_path, 'r') as f:
            unsolved_datas = filter(lambda arg: arg[0] is not None,
                                    [DataDownloadTool.format_data(line) for line in
                                     filter(lambda arg: arg.startswith('>'), f.readlines())])
            unsolved_datas = list(unsolved_datas)
        for iteration_time in range(1, ExperimentConfig.MAX_ITERATION_TIME + 1):
            unsolved_datas = self.find_neighbor_batch(unsolved_datas, iteration_time)
            if len(unsolved_datas) == 0:
                break
        print("Failed data:" + str(len(unsolved_datas)) + "," + str(unsolved_datas))
        self.source_gene_distribution_analysis()

    @staticmethod
    def analysis_download_file(download_file_path, inter):
        left = min(inter)
        right = max(inter)
        gene_info = NCBIDatabase(download_file_path)
        if not gene_info.initialize():
            return False, None
        near_small = None
        near_big = None
        res_set = set()
        for idx, gene_segment in enumerate(gene_info.gene_segments):
            if gene_segment.cds[1] <= left:
                if not near_small or near_small.cds[1] < gene_segment.cds[1]:
                    near_small = gene_segment
            if gene_segment.cds[0] >= right:
                if not near_big or near_big.cds[0] > gene_segment.cds[0]:
                    near_big = gene_segment
            if gene_segment.cds[0] <= left <= gene_segment.cds[1]:
                res_set.add(str(gene_segment))
            if gene_segment.cds[0] <= right <= gene_segment.cds[1]:
                res_set.add(str(gene_segment))
        if near_small:
            res_set.add(near_small)
        if near_big:
            res_set.add(near_big)
        sequence = gene_info.dna_code[left - 1:right]
        if inter[0] > inter[1]:
            sequence = get_opposite_dna(sequence[::-1])
        return True, {'source': gene_info.source, 'data': list(res_set), 'sequence': sequence}

    @staticmethod
    def check_gene(left, right, direction, gene, target):
        if direction == '-':
            peer = min(left, right)
            gene_peer = max(gene['left'], gene['right'])
            if peer > gene_peer:
                return target is None or max(target['left'], target['right']) < gene_peer
        elif direction == '+':
            peer = max(left, right)
            gene_peer = min(gene['left'], gene['right'])
            if peer < gene_peer:
                return target is None or min(target['left'], target['right']) > gene_peer
        else:
            raise ValueError('Direction should be - or +')

    @staticmethod
    def read_gene(line):
        inter, gene = line.strip().split('\t')
        left, right = map(int, inter.split('-'))
        return {
            'gene': gene,
            'left': left,
            'right': right
        }
