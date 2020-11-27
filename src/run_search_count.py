import os
from urllib import request

from experiment_config import ExperimentConfig
import re

input_file_name = '20201126 0vs05_NCBI_search.txt'
url = 'https://pubmed.ncbi.nlm.nih.gov/?term='


def extract_result_info(file_path):
    start = '<div class="results-amount">'
    step = -1
    for line in open(file_path, 'rb'):
        line = line.decode('utf8', errors='ignore')
        line = line.strip()
        if line == start:
            step = 0
        elif step == 0 and line:
            if line.find('No results were found') >= 0:
                return '0'
            else:
                return line.replace('<span class="value">', '').replace('</span>', '')
    return '1'


def check_file(file_path):
    if not os.path.exists(file_path):
        return False
    with open(file_path, 'rb') as fr:
        buff = fr.read().strip()
        return len(buff) > 0


def download_gene_key(gene, key, file_path):
    if check_file(file_path):
        return True
    try:
        x = request.urlopen(url + gene + '+' + key, timeout=30)
        body = x.read()
        with open(file_path, 'wb') as fw:
            fw.write(body)
        return True
    except Exception as e:
        print(e)
        return False


def download_gene(gene):
    result = [gene]
    for key in ['drug', 'Kanamycin+B', 'aminoglycoside', 'antibiotic', 'biofilm']:
        key = key.replace('+', '_')
        file_path = os.path.join(ExperimentConfig.pubmed_ncbi_directory,
                                 '{0}_{1}.html'.format(gene, key))
        retry = 0
        while not download_gene_key(gene, key, file_path) and retry < 3:
            retry += 1
        if retry == 3:
            print("Error: {0}_{1}".format(gene, key))
            result.append('Error')
        else:
            print("Success: {0}_{1}".format(gene, key))
            result.append(extract_result_info(file_path))
    return result


def download_all(gene_list):
    with open(os.path.join(ExperimentConfig.output_directory, 'pubmed_count.tsv'), 'w') as fw:
        fw.write('\t'.join(['gene', 'drug', 'Kanamycin+B', 'aminoglycoside', 'antibiotic', 'biofilm']) + '\n')
        for gene in gene_list:
            result = download_gene(gene.strip())
            fw.write('\t'.join(result) + '\n')
            fw.flush()


if __name__ == '__main__':
    with open(os.path.join(ExperimentConfig.data_directory, input_file_name), 'r') as fr:
        gene_list = fr.readlines()[1:]
        download_all(gene_list)
