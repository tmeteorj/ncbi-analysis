import os
import unittest
from urllib import request

from analysis.cluster_match import ClusterMatcher
from analysis.gene_extract import GeneExtract
from analysis.neighbor_analysis import NeighborAnalysis
from experiment_config import ExperimentConfig
from utils.data_download_util import DataDownloadTool
from utils.gene_file_util import GeneSegment, GeneFileReader


class TestRNAAnalysis(unittest.TestCase):
    def setUp(self):
        self.root_directory = os.sep.join(os.getcwd().split(os.sep)[:-2])
        self.data_directory = os.path.join(self.root_directory, 'data', 'rna_analysis')
        self.output_directory = os.path.join(self.root_directory, 'data', 'rna_analysis_result')
        self.download_directory = os.path.join(self.data_directory, 'rna_download_data')
        self.test_fna_name = '18_utr_bacteria.fna'
        self.rna_tag = '18_utr'

    def test_cluster_matcher(self):
        input_path = os.path.join(self.data_directory, self.test_fna_name)
        cluster_match = ClusterMatcher(self.rna_tag, input_path, self.output_directory)
        cluster_match.run()

    def test_neighbor_analysis(self):
        file_name = self.rna_tag + '_all_result.txt'
        input_path = os.path.join(self.output_directory, file_name)
        neighbor_analysis = NeighborAnalysis(input_path, self.download_directory, self.output_directory)
        neighbor_analysis.run()

    def test_gene_segment(self):
        gene_segment = GeneSegment()
        gene_segment.extract_attribute('/product=\"asd')
        self.assertEqual(gene_segment.product, 'asd')
        gene_segment.extract_attribute('/product=\"asdasd\"')
        self.assertEqual(gene_segment.product, 'asdasd')

    def test_download_tool(self):
        key = "NZ_FTWJ01000011.1"
        url = ExperimentConfig.URL_LIB_PREFIX + key
        res = request.urlopen(url, timeout=60)
        self.assertIsNotNone(res)
        output_file_path = os.path.join(self.output_directory, 'sample_download_%s.txt' % key)
        self.assertTrue(DataDownloadTool.download_data(key, output_file_path))

    def test_gene_data_reader(self):
        input_path = os.path.join(self.download_directory, 'NC_000913.3.txt')
        gene_data_reader = GeneFileReader(input_path)
        gene_data_reader.build_information()
        self.assertTrue(len(gene_data_reader.gene_segments) > 0)
        with open(os.path.join(self.data_directory, 'gene_all.txt'), 'w', encoding='utf8') as fw:
            for gene_segment in gene_data_reader.gene_segments:
                if gene_segment.gene is not None:
                    fw.write(gene_segment.gene + '\n')

    def test_extract_gene(self):
        data_path = os.path.join(self.download_directory, 'NC_000913.3.txt')
        rna_path = os.path.join(self.data_directory, 'rna_sequence.txt')
        output_path = os.path.join(self.output_directory, 'rna_sequence_result.tsv')
        tool = GeneExtract(data_path, rna_path, output_path)
        tool.run()
