from analysis.ecocyc_analysis import EcocycAnalysis
from analysis.gene_extract import GeneExtract
from experiment_config import *

extract_gene_file_name = 'NC_000913.3.txt'
ecocyc_gene_files = ['20210609.txt']
ecocyc_params = {
    'from_gene_names': True,
    'output_best_promoter': True,
    'output_gene_sequence': False,
    'output_detail_information': True,
    'analysis_promoter': True,
    'if_get_summary': True,
    'if_get_go_table': True
}
from_gene_names = True
output_best_promoter = True
cookie = 'JSESSIONID=3D7EA386ED826530EA8854B46360216E; _ga=GA1.2.1085079220.1623245828; _gid=GA1.2.442820656.1623245828; PTools-session=biocyc13a~biocyc13-3832234589|NIL NIL NIL NIL 0 (:WEB NIL -1 ((:BASICS -1) (:QUERIES -1) (:ADVANCED -1))) NIL NIL NIL ECOBASE NIL NIL |6a977duko82vj71z42jcw7y5vjmmi0d; frameHeight=695; frameWidth=623; pagecount=8; _gat=1; _gat_UA-7250911-5=1'
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
