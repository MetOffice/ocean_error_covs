"""
Top script to compute Desroziers statistics for a given month

@author: A Fowler Nov 2025

Modified by Y Chen, University of Reading, 2025
"""
import sys

from .config import Config
from .estimator import CovStatBinner, CorrEstimator
from . import log
from .grid import Grid
from .obs_factory import ObsFactory
from .predictor import Predictor
from .output_handler import OutputHandler

def calc_desroziers(configfile='config.ini', mask_func=None):
    """Calculate uncertainties from Desroziers statistics.

    Information should be given by configuration files.

    Parameters
    ----------
    configfile : str, optional
        Path to the top-level configuration file. Defaults to ``'config.ini'``
        in the current working directory.
    mask_func : callable, optional
        A function that selects observations belonging to a given predictor bin.
        If not provided, observations are selected by matching the ``predictor``
        variable in the input data to the bin index (or all observations are
        selected when no predictor variable is present).

        The function must have the signature::

            def mask_func(name: str, data: dict[str, np.ndarray], i_bin: int) -> np.ndarray

        Parameters
        ~~~~~~~~~~
        name : str
            Name of the observation type as specified in the ``[General]``
            section of the observation ``.ini`` file.
        data : dict[str, np.ndarray]
            Dictionary of all variables loaded for the current time step.
            Always contains at least ``'lon'`` and ``'lat'``.  Any additional
            variables listed under ``variables`` in the ``.ini`` file are also
            present, including the predictor variable specified by ``pred_name``.
        i_bin : int
            Zero-based index of the predictor bin to select.

        Returns
        ~~~~~~~
        np.ndarray
            Boolean array of length ``len(data['lon'])`` where ``True``
            indicates observations that belong to bin ``i_bin``.

        Example
        ~~~~~~~
        The following example selects observations by dominant Optical Water
        Type (OWT), while excluding observations where the total OWT membership
        score is outside a quality range::

            import numpy as np

            def mask_func_owt(name, data, i_bin):
                total_owt_limits = [0.5, 1.5]
                mask = data['OWT_dom'].copy()
                mask[data['OWT_tot'] < total_owt_limits[0]] = -1
                mask[data['OWT_tot'] > total_owt_limits[1]] = -1
                return mask == i_bin

            calc_desroziers(mask_func=mask_func_owt)
    """
    # Load configuration
    config_t = Config(configfile=configfile)
    config_t.load_obs_configs()

    estimates = config_t['Uncertainty'].get('estimate', 'R,HBH,R+HBH').split(',')
    is_horizontal = config_t['Uncertainty'].getboolean('is_horizontal', True)

    # create an observation object
    obs_factory = ObsFactory(config_t.obs_configs, estimates)
    # initialise output handler
    output = OutputHandler(config_t['Output'])
    # Create components
    for obs_configs in obs_factory.obs_pair_iterator(config_t.obs_configs):
        list_obs_types = list(obs_configs.keys())
        predictor = Predictor(obs_configs, mask_func=mask_func)
        grid = Grid(obs_configs[list_obs_types[0]]['Grid'], is_horizontal)
        # Create estimator with only what it needs
        estimator = CovStatBinner(grid, predictor, list_obs_types,
                                  estimates)
        # Pass data iterator to estimator
        results = estimator.do_binning(obs_factory.data_iterator(list_obs_types, is_horizontal))
        # Save output
        output.save_binned(predictor, grid, results)

        # Post-process results
        ss_lim = obs_configs[list_obs_types[0]]['Predictor'].getint('n_samples', 100)
        smooth_window = obs_configs[list_obs_types[0]]['Predictor'].getint('smooth_window', 6)
        cov_estimator = CorrEstimator(estimates, ss_lim=ss_lim, smooth_window=smooth_window)
        cov_estimates = cov_estimator.do_compute(results)
        output.save_cov(predictor, grid, cov_estimates, list_obs_types)

    log.logger.info('*************************************')
    msg = sys.argv[0]+'*** COMPLETED ***'
    log.logger.info(msg)
    log.logger.info('*************************************')


def _cli():
    """Console-script entry point for the ``calc_desroziers`` command."""
    configfile = sys.argv[1] if len(sys.argv) > 1 else 'config.ini'
    calc_desroziers(configfile=configfile)

if __name__ == "__main__":
    _configfile = sys.argv[1] if len(sys.argv) > 1 else 'config.ini'
    calc_desroziers(configfile=_configfile)