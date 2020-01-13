import os
import re
import sys

from experiment_config import ExperimentConfig
from utils.factories.logger_factory import LoggerFactory
from utils.str_util import StrConverter


class ClusterMatcher:
    def __init__(self, rna_tag, input_path, output_directory):
        self.rna_tag = rna_tag
        self.input_path = input_path
        self.output_directory = output_directory
        self.logger = LoggerFactory(1)

        file_prefix = StrConverter.extract_file_name(rna_tag)
        self.cluster_result_path = os.path.join(self.output_directory,
                                                '%s_cluster_result.txt' % file_prefix)
        self.sample_result_path = os.path.join(self.output_directory,
                                               '%s_sample_result.txt' % file_prefix)
        self.all_result_path = os.path.join(self.output_directory,
                                            '%s_all_result.txt' % file_prefix)
        self.only_result_path = os.path.join(self.output_directory,
                                             '%s_only_result.txt' % file_prefix)

    def run(self):
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        data = self.read_data()
        self.logger.info('data size = %d' % len(data))
        self.compare_data(data)

    def format_data(self, index, lines):
        ig, gene_no = self.should_ignore(index, lines[3])
        if ig and not gene_no:
            return None
        dic = {'ignored': ig, 'index': index, 'geneNo': gene_no}
        if dic['ignored']:
            return dic

        data = [{}, {}, {}]
        action = 0
        others = lines[:4]
        for line in lines[4:]:
            if line.strip() == '':
                continue
            if line.strip().startswith(self.rna_tag):
                action = 1
                self.update_sequence(index, data[0], line)
            elif action == 1:
                action = 2
                self.update_sequence(index, data[1], line)
            elif action == 2:
                action = 0
                self.update_sequence(index, data[2], line)
            else:
                action = 0
        dic['data'] = data
        dic['others'] = others
        return dic

    def read_data(self):
        data = []
        buff = []
        index = 0
        for line in open(self.input_path, 'r'):
            if line.startswith('>>'):
                if len(buff) > 0:
                    index += 1
                    dic = self.format_data(index, buff)
                    if dic:
                        data.append(dic)
                buff = []
            buff.append(line)
        if len(buff) > 0:
            index += 1
            data.append(self.format_data(index, buff))
        return data

    def compare_data(self, data):
        same = {}
        cluster = {}
        si = len(data)
        sample_cluster = {}
        for i in range(si):
            if i in same or data[i]['ignored']:
                continue
            same[i] = i
            cluster[i] = [data[i]['geneNo']]
            sample_cluster[i] = [data[i]]
            for j in range(i + 1, si):
                if j in same or data[j]['ignored']:
                    continue
                if data[i]['data'][1]['seq'].upper() == data[j]['data'][1]['seq'].upper():
                    same[j] = i
                    cluster[i].append(data[j]['geneNo'])
                    sample_cluster[i].append(data[j])
        fw = open(self.cluster_result_path, 'w')
        for _ in cluster:
            fw.write('%d\t%s\n' % (len(cluster[_]), ','.join(map(str, cluster[_]))))
        fw.close()
        fw = open(self.sample_result_path, 'w')
        for _ in sample_cluster:
            sample = sample_cluster[_][0]
            fw.write(''.join(sample['others']))
            fw.write('\n')
            for elem in sample['data']:
                fw.write('%19s %8s %131s %8s\n' % (
                    elem.get('name', ''), elem.get('start', ''), elem.get('seq', ''), elem.get('end', '')))
            fw.write('\n')
        fw.close()
        fw_all = open(self.all_result_path, 'w')
        fw_only = open(self.only_result_path, 'w')
        other = set()
        for _ in sample_cluster:
            for item in sample_cluster[_]:
                elem = item['data'][-1]
                flag = True
                for x in elem['seq'].strip():
                    if x.upper() in set('AUCG'): continue
                    other.add(x.upper())
                    flag = False
                fw_all.write('>%s/%s-%s\n%s\n' % (elem['name'], elem['start'], elem['end'], elem['seq'].upper()))
                fw_all.write('\n')
                if not flag:
                    continue
                fw_only.write('>%s/%s-%s\n%s\n' % (elem['name'], elem['start'], elem['end'], elem['seq'].upper()))
                fw_only.write('\n')
            fw_all.write('\n')
            fw_only.write('\n')
        fw_all.close()
        fw_only.close()
        self.logger.info('\n'.join(list(other)))

    def should_ignore(self, index, line):
        info = re.split(r'\s+', line.strip())
        gene_no = info[0].strip('()')
        if info[1] == '?':
            return False, gene_no
        elif info[1] == '!':
            return False, gene_no
        else:
            self.logger.info('ignore check failed: %d' % index)
            return True, None

    def update_sequence(self, index, elem, line):
        if line.strip()[-1] not in ExperimentConfig.SET_NUMBER_RANGE10:
            elem['seq'] = elem.get('seq', '') + line.strip()
            return
        try:
            info = re.split(r'\s+', line.strip())
            name = info[0]
            start = int(info[1])
            end = int(info[-1])
            seq = ' '.join(info[2:-1])
        except:
            self.logger.info('value num is not ok: %d, %s' % (index, line))
            sys.exit(1)
        if elem.get('name', name) != name:
            self.logger.info('name is not equal: %d' % index)
            sys.exit(1)
        elem['name'] = name
        start = int(start)
        end = int(end)
        if 'start' not in elem:
            elem['start'] = start
        elem['end'] = end
        elem['seq'] = elem.get('seq', '') + seq
