import os

from analysis.cluster_match import ClusterMatcher
from analysis.gene_extract import GeneExtract
from analysis.gene_name_eco_download import EcocycAnalysis
from analysis.neighbor_analysis import NeighborAnalysis

root_directory = os.sep.join(os.getcwd().split(os.sep)[:-1])
data_directory = os.path.join(root_directory, 'data', 'rna_analysis')
output_directory = os.path.join(root_directory, 'data', 'rna_analysis_result')
rna_download_directory = os.path.join(data_directory, 'rna_download_data')
ecocyc_download_directory = os.path.join(data_directory, 'ecocyc_download_data')

fna_name = 'T2_sun_bacteria.fna'
rna_tag = 'T2_sun'
do_cluster_match = False
do_neighbor_analysis = False

extract_gene_file_name = 'NC_000913.3.txt'
extract_gene_sequence = 'rna_sequence.txt'
do_gene_extract = False

ecocyc_gene_files = ['transcription_ul.txt', 'tRNA_ul.txt']
from_gene_names = False
do_ecocyc_analysis = True


def run_cluster_match():
    input_path = os.path.join(data_directory, fna_name)
    cluster_match = ClusterMatcher(rna_tag, input_path, output_directory)
    cluster_match.run()


def run_neighbor_analysis():
    file_name = rna_tag + '_all_result.txt'
    input_path = os.path.join(output_directory, file_name)
    neighbor_analysis = NeighborAnalysis(input_path, rna_download_directory, output_directory)
    neighbor_analysis.run()
    # neighbor_analysis.source_gene_distribution_analysis()


def run_gene_extract():
    data_path = os.path.join(rna_download_directory, extract_gene_file_name)
    rna_path = os.path.join(data_directory, 'rna_sequence.txt')
    output_path = os.path.join(output_directory, 'rna_sequence_result.tsv')
    gene_extract = GeneExtract(data_path, rna_path, output_path)
    gene_extract.run()


def run_ecocyc_analysis():
    for ecocyc_gene_file in ecocyc_gene_files:
        input_path = os.path.join(data_directory, ecocyc_gene_file)
        ecocyc_analysis = EcocycAnalysis(input_path, ecocyc_download_directory, output_directory, from_gene_names)
        ecocyc_analysis.run()


if __name__ == '__main__':
    if do_cluster_match:
        run_cluster_match()
    if do_neighbor_analysis:
        run_neighbor_analysis()
    if do_gene_extract:
        run_gene_extract()
    if do_ecocyc_analysis:
        run_ecocyc_analysis()
