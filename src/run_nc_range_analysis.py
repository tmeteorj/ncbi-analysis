from analysis.gene_range_analysis import GeneRangeExtract
from experiment_config import *

nc_gene_file_name = 'NC_000913.3.txt'

if __name__ == '__main__':
    data_path = os.path.join(ExperimentConfig.rna_download_directory, nc_gene_file_name)
    gene_range_extract = GeneRangeExtract(data_path, ExperimentConfig.output_directory)
    gene_range_extract.run()
