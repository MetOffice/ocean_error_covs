"""General configuration for Desroziers covariance estimation

Author: Y Chen, University of Reading, 2025
"""

import configparser
import os

from . import log

class Config:
    """This class controls the configuration/namelist reading
    """
    _initialised = False
    def __init__(self, configfile='config.ini'):
        assert not Config._initialised, 'Configuration is already set up!!!'
        Config._initialised = True

        self.config = configparser.ConfigParser()
        self.config.read(configfile)
        self._configfile = os.path.abspath(configfile)

    def load_obs_configs(self) -> None:
        obs_config_fnames = self.config['ObsTypes'
                                        ].get('obs_config_files',''
                                              ).replace(' ','').split(',')
        config_dir = os.path.dirname(self._configfile)
        self.obs_configs = {}
        for obs_config_fname in obs_config_fnames:
            candidates = [
                obs_config_fname,
                obs_config_fname + '.ini',
                os.path.join(config_dir, obs_config_fname),
                os.path.join(config_dir, obs_config_fname + '.ini'),
                os.path.abspath(obs_config_fname),
                os.path.abspath(obs_config_fname + '.ini'),
            ]
            path = next((c for c in candidates if os.path.isfile(c)), None)
            if path is None:
                raise FileNotFoundError(
                    f"Observation config file '{obs_config_fname}' not found.\n"
                    f"Searched (with and without '.ini' suffix):\n"
                    f"  {config_dir} (directory of config file)\n"
                    f"  {os.getcwd()} (current working directory)\n"
                    "Provide an absolute path in the config file."
                )
            log.logger.info(f"Loading observation config from: {path}")
            parser = configparser.ConfigParser()
            parser.read(path)
            self.obs_configs[obs_config_fname] = parser

    def __getitem__(self, key):
        return self.config[key]
