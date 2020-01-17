from analysis.gene_extract import GeneExtract
from analysis.gene_name_eco_download import EcocycAnalysis
from experiment_config import *

extract_gene_file_name = 'NC_000913.3.txt'
ecocyc_gene_files = ['77.txt']
ecocyc_params = {
    'from_gene_names': True,
    'output_best_promoter': True,
    'output_gene_sequence': True,
    'output_detail_information': True
}
from_gene_names = True
output_best_promoter = True
cookie = 'pagecount=6; JSESSIONID=C53ACDDDFED95199673F71FD85C1E3D3; _gat=1; _gid=GA1.2.1786139794.1579272089; PTools-session=biocyc14b~biocyc14-3786098971%7CNIL%20NIL%20%22%22%20NIL%200%20(%3AWEB%20NIL%20-1%20((%3ABASICS%20-1)%20(%3AQUERIES%20-1)%20(%3AADVANCED%20-1)))%20NIL%20NIL%20ECOBASE%20NIL%20NIL%20%7Ch0x9lugmbx6fhk3bevsf98cen2omar5; _ga=GA1.2.407871027.1577110083; windowOrg=ptools0%3AECOLI%3A; recentOrgID0=ECOLI; frameWidth=1500; frameHeight=761'

if __name__ == '__main__':
    data_path = os.path.join(ExperimentConfig.rna_download_directory, extract_gene_file_name)
    for ecocyc_gene_file in ecocyc_gene_files:
        input_path = os.path.join(ExperimentConfig.data_directory, ecocyc_gene_file)
        ecocyc_analysis = EcocycAnalysis(input_path, ExperimentConfig.ecocyc_download_directory,
                                         ExperimentConfig.output_directory, ecocyc_params, cookie)
        ecocyc_analysis.run()
        if ecocyc_params['output_gene_sequence']:
            seq_path = ecocyc_analysis.ecocyc_result_path
            left_idx = ecocyc_analysis.sequence_start_idx
            right_idx = ecocyc_analysis.sequence_end_idx
            gene_extract = GeneExtract(data_path, seq_path, ExperimentConfig.output_directory, 'range', left_idx,
                                       right_idx)
            gene_extract.run()
