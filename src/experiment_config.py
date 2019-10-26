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
VALUE_REPEAT_REGION_START = 'repeat_region '
VALUE_GENE_START_WITH_SPACE = '     gene'
VALUE_CDS_START_WITH_SPACE = '     CDS'

SET_NUMBER_RANGE10 = set(map(str, range(10)))

MAX_ITERATION_TIME = 100
