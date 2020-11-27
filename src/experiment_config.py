import os


class ExperimentConfig:
    URL_LIB_PREFIX = 'https://www.ncbi.nlm.nih.gov/nuccore/'
    URL_DATA_DOWNLOAD = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=%s&db=nuccore' \
                        '&report=genbank&extrafeat=976&conwithfeat=on&withparts=on&retmode=txt' \
                        '&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=%s'
    MAX_DOWNLOAD_RETRY_TIME = 3

    KEY_LIB_UID = 'ncbi_uidlist'
    KEY_GENE = 'gene'
    KEY_PRODUCT = 'product'
    VALUE_UNKNOWN = 'UNKNOWN'

    VALUE_SOURCE_START = 'SOURCE'
    VALUE_GENE_START = 'gene    '
    VALUE_CODE_START = 'CDS    '
    VALUE_DNA_PART_END = '//'
    # Ignore repeat region
    VALUE_REPEAT_REGION_START = 'repeat_region '
    # Not ignore repeat region
    VALUE_GENE_START_WITH_SPACE = '     gene'
    VALUE_CDS_START_WITH_SPACE = '     CDS'

    SET_NUMBER_RANGE10 = set(map(str, range(10)))

    MAX_ITERATION_TIME = 100
    MAX_THREAD_NUM = 16

    root_directory = os.sep.join(os.getcwd().split(os.sep)[:-1])
    data_directory = os.path.join(root_directory, 'data', 'rna_analysis')
    output_directory = os.path.join(root_directory, 'data', 'rna_analysis_result')
    rna_download_directory = os.path.join(data_directory, 'rna_download_data')
    ecocyc_download_directory = os.path.join(data_directory, 'ecocyc_download_data')
    pubmed_ncbi_directory = os.path.join(data_directory, 'pubmed_ncbi')

    @staticmethod
    def disable_repeat_region():
        ExperimentConfig.VALUE_REPEAT_REGION_START = 'DEPRECATED: repeat_region '

    @staticmethod
    def enable_repeat_region():
        ExperimentConfig.VALUE_REPEAT_REGION_START = 'repeat_region '
