import os

from analysis.cluster_match import ClusterMatcher
from analysis.gene_extract import GeneExtract
from analysis.gene_name_eco_download import EcocycAnalysis
from analysis.gene_range_analysis import GeneRangeExtract
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
extract_gene_sequence = 'gene_utr_to_find.txt'
gene_extract_based = 'range'
do_gene_extract = False

ecocyc_gene_files = ['gene_all.txt']
ecocyc_params = {
    'from_gene_names': True,
    'output_best_promoter': True,
    'output_gene_sequence': True,
    'output_detail_information': True
}
from_gene_names = True
output_best_promoter = True
cookie = 'windowOrg=ptools0%3AECOLI%3Aptools1%3AECOLI%3A; recentOrgID0=ECOLI; pagecount=16; frameHeight=759; JSESSIONID=1D257C680D66CA41ED0FC9432FD0BD3E; frameWidth=1500; _gat=1; _gid=GA1.2.610587971.1577521867; PTools-session=biocyc14b~biocyc14-3786098971%7CNIL%20NIL%20%22%22%2042107%200%20(%3AWEB%20NIL%20-1%20((%3ABASICS%20-1)%20(%3AQUERIES%20-1)%20(%3AADVANCED%20-1)))%20NIL%20NIL%20ECOBASE%20NIL%20NIL%20%7Cd0vt8w4rc07g29ilyr3o518e74ero6c; _ga=GA1.2.407871027.1577110083; credentialId=218865; secretKey=27oXO8IVRHh01SA3ae/qL9Yqfwk='
do_ecocyc_analysis = False

nc_gene_file_name = 'NC_000913.3.txt'
do_nc_range_extract = True


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
    rna_path = os.path.join(data_directory, extract_gene_sequence)
    gene_extract = GeneExtract(data_path, rna_path, output_directory, gene_extract_based)
    gene_extract.run()


def run_ecocyc_analysis():
    data_path = os.path.join(rna_download_directory, extract_gene_file_name)
    for ecocyc_gene_file in ecocyc_gene_files:
        input_path = os.path.join(data_directory, ecocyc_gene_file)
        ecocyc_analysis = EcocycAnalysis(input_path, ecocyc_download_directory, output_directory, ecocyc_params, cookie)
        ecocyc_analysis.run()
        if ecocyc_params['output_gene_sequence']:
            seq_path = ecocyc_analysis.ecocyc_result_path
            left_idx = ecocyc_analysis.sequence_start_idx
            right_idx = ecocyc_analysis.sequence_end_idx
            gene_extract = GeneExtract(data_path, seq_path, output_directory, 'range', left_idx, right_idx)
            gene_extract.run()


def run_nc_range_extract():
    data_path = os.path.join(rna_download_directory, nc_gene_file_name)
    gene_range_extract = GeneRangeExtract(data_path, output_directory)
    gene_range_extract.run()


if __name__ == '__main__':
    if do_cluster_match:
        run_cluster_match()
    if do_neighbor_analysis:
        run_neighbor_analysis()
    if do_gene_extract:
        run_gene_extract()
    if do_ecocyc_analysis:
        run_ecocyc_analysis()
    if do_nc_range_extract:
        run_nc_range_extract()
