"""Handling different observation types

Author: Y Chen, University of Reading, 2025
"""
from collections import OrderedDict
import configparser
import itertools
import typing

import numpy as np

from .input_handler import InputHandler
from . import log

class ObsFactory:
    """Class to handle different observation types based on configuration.
    """
    def __init__(self, configs:dict[str, configparser.ConfigParser],
                 estimates: list[str]) -> None:
        log.logger.info('Setting up observation types...')

        self.names = []
        self.long_names = []
        self.inputs = OrderedDict()
        self.obs_name_cfg_map = OrderedDict()
        for obs_config_name, config in configs.items():
            self.names.append(
                config['General'].get('name', obs_config_name)
            )
            self.obs_name_cfg_map[self.names[-1]] = obs_config_name
            self.long_names.append(
                config['General'].get('long_name', self.names[-1])
            )
            self.inputs[self.names[-1]] = InputHandler(config['Input'], estimates)

        # get pairs of observation types for covariance calculations
        self.obs_pairs = itertools.product(self.names, repeat=2)


    def obs_pair_iterator(self, configs:dict[str, configparser.ConfigParser]
                          ) -> typing.Iterator[dict[str, configparser.ConfigParser]]:
        """Iterator over observation type pairs.

        Yields
        ------
        dict[str, configparser.ConfigParser]
            A dictionary containing a pair of observation type configurations.
        """
        for obs_a, obs_b in self.obs_pairs:
            if obs_a == obs_b:
                msg = f'Processing covariance for observation type: {obs_a}'
            else:
                msg = f'Processing cross-covariance for observation types: {obs_a}, {obs_b}'
            log.logger.info(msg)
            obs_configs = OrderedDict()
            obs_configs[obs_a] = configs[self.obs_name_cfg_map[obs_a]]
            obs_configs[obs_b] = configs[self.obs_name_cfg_map[obs_b]]
            yield obs_configs

    def data_iterator(self, names:list[str], is_horizontal: bool
                      ) -> typing.Iterator[dict[str, dict[str, np.ndarray]]]:
        """Iterator over NetCDF files in the specified mode.

        Yields
        ------
        dict[str, dict[str, np.ndarray]]
            A dictionary of variable names and their corresponding data arrays.
        """
        # iterator for each data file.
        self._data_iterator = [self.inputs[name].read(is_horizontal) for name in names]

        for v_dict in zip(*self._data_iterator):
            data:dict[str, dict[str, np.ndarray]] = OrderedDict()
            for name, d in zip(names, v_dict):
                data[name] = d
            yield data