import os
import unittest

from experiment_config import ExperimentConfig
from utils.atcc_database import ATCCDatabase


class TestColiDatabase(unittest.TestCase):
    def test_load_database(self):
        file_path = os.path.join(ExperimentConfig.rna_download_directory, 'CP009072.1.txt')
        coli_database = ATCCDatabase(file_path)
        self.assertTrue(len(coli_database.gene_segments) > 0)
