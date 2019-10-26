import os

from analysis.cluster_match import ClusterMatcher
from analysis.gene_extract import GeneExtract
from analysis.neighbor_analysis import NeighborAnalysis

root_directory = os.sep.join(os.getcwd().split(os.sep)[:-1])
data_directory = os.path.join(root_directory, 'data', 'rna_analysis')
output_directory = os.path.join(root_directory, 'data', 'rna_analysis_result')
download_directory = os.path.join(data_directory, 'rna_download_data')

fna_name = '7a_utr_bacteria.fna'
rna_tag = '7a_utr'
extract_gene_file_name = 'NC_000913.3.txt'
extract_gene_sequence = 'rna_sequence.txt'
do_cluster_match = True
do_neighbor_analysis = True
do_gene_extract = True


def run_cluster_match():
    input_path = os.path.join(data_directory, fna_name)
    cluster_match = ClusterMatcher(rna_tag, input_path, output_directory)
    cluster_match.run()


def run_neighbor_analysis():
    file_name = rna_tag + '_all_result.txt'
    input_path = os.path.join(output_directory, file_name)
    neighbor_analysis = NeighborAnalysis(input_path, download_directory, output_directory)
    neighbor_analysis.run()


def run_gene_extract():
    data_path = os.path.join(download_directory, extract_gene_file_name)
    rna_path = os.path.join(data_directory, 'rna_sequence.txt')
    output_path = os.path.join(output_directory, 'rna_sequence_result.tsv')
    gene_extract = GeneExtract(data_path, rna_path, output_path)
    gene_extract.run()


if __name__ == '__main__':
    if do_cluster_match:
        run_cluster_match()
    if do_neighbor_analysis:
        run_neighbor_analysis()
    if do_gene_extract:
        run_gene_extract()
