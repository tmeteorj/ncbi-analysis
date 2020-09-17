import os
from abc import ABC

from experiment_config import ExperimentConfig


class BaseExperiment(ABC):
    input_file: str
    experiment_name: str = 'base_experiment'
    working_directory: str = None

    def __init__(self, input_file: str):
        self.input_file = input_file

