import os

from analysis.found_gene_name_near_a_range import FoundGeneNameNearARange
from experiment_config import ExperimentConfig
from utils.gene_position_helper import GenePositionHelper
from utils.ncbi_database import NCBIDatabase

input_file_name = 'find_gene\\input_with_locus.txt'
data_name = 'NC_000913.3.txt'

if __name__ == '__main__':
    database = NCBIDatabase(os.path.join(ExperimentConfig.rna_download_directory, data_name))
    helper = GenePositionHelper(database)
    exp = FoundGeneNameNearARange(helper)
    exp.run(os.path.join(ExperimentConfig.data_directory, input_file_name))
