"""
Output handler for saving Desroziers statistics.

Author: Y Chen, University of Reading, 2025
"""
import configparser
import os
import time
import re

import netCDF4
import numpy as np

from .predictor import Predictor
from .grid import Grid
from .estimator import BinnedStats, CovarianceEstimates

class OutputHandler:
    def __init__(self, config_t:configparser.SectionProxy):
        self.outputdir = config_t.get('output_dir','./')
        os.makedirs(self.outputdir, exist_ok=True)
        self.binned_fname = config_t.get('binned_filename','binned_data.nc')
        self.cov_fname = config_t.get('cov_filename','cov_data.nc')
        self.binned_path = os.path.join(self.outputdir,self.binned_fname)
        self.cov_path = os.path.join(self.outputdir,self.cov_fname)

    def save_binned(self, predictor:Predictor, grid:Grid, results:BinnedStats) -> None:
        """
        Save binned statistics to a NetCDF file.

        Parameters
        ----------
        predictor : Predictor
            The predictor object containing binning information.
        grid : Grid
            The grid object containing binning information.
        results : BinnedStats
            The binned statistics to be saved.

        Returns
        -------
        None
        """
        # write data to file
        if len(results.n_sample_obs) == 1:
            obs_a, = list(results.n_sample_obs.keys())
            obs_b = obs_a
        else:
            (obs_a, obs_b) = results.n_sample_obs.keys()
        fname = f'{self.binned_path}_{obs_a}-{obs_b}.nc'
        dataset = netCDF4.Dataset(fname,'w',format='NETCDF4_CLASSIC') # pylint: disable=no-member
        # Global Attributes
        dataset.description = f'Dataset containing outputted innovation and residual statistics binned by {predictor.name}'
        dataset.history = f'Created {time.ctime(time.time())}'

        dataset.createDimension('n_grid', grid.n_gp)
        dataset.createDimension('n_predictors', predictor.n_classes)

        total_obs = dataset.createVariable('total_obs_bin',
                                            np.int32,
                                            ('n_predictors','n_grid')
                                            )
        total_obs.setncatts({'long_name': f"total number of observations ({obs_a}x{obs_b})"
                             " in each predictor bin"})
        total_obs[:] = results.n_obs_bin

        bins = dataset.createVariable('bins',np.float32,('n_grid',))
        bins.setncatts({'long_name': "Grid of correlation length", 'units': "km"})
        bins[:] = grid.d


        longnames = {'d_ob_i': "Sum of y-H(xb) for starting gridpoint of spatial correlation",
                     'd_ob_j': "Sum of y-H(xb) for end gridpoint of spatial correlation",
                     'd_oa': "Sum of y-H(xa) for end gridpoint of spatial correlation",
                     'd_ab': "Sum of H(xa)-H(xb) for end gridpoint of spatial correlation",
        }
        for v in results.d_sum:
            nc_var = dataset.createVariable(v, np.float32, ('n_predictors','n_grid'))
            nc_var.setncatts({'long_name': longnames[v], 'units': "none"})
            nc_var[:] = results.d_sum[v]

        longnames = {'HBH': "Sum of (y-H(xb))*(H(xa)-H(xb).T",
                     'R': "Sum of (y-H(xb))*(y-H(xa)).T",
                     'R_HBH': "Sum of (y-H(xb))*(y-H(xb)).T",
                     }
        for v in results.err_sum:
            nc_var = dataset.createVariable(v, np.float32, ('n_predictors','n_grid'))
            nc_var.setncatts({'long_name': longnames[v], 'units': "none"})
            nc_var[:] = results.err_sum[v]

        dataset.close()

    def save_cov(self, predictor:Predictor, grid:Grid, results:CovarianceEstimates,
                 obs_names:list[str]) -> None:
        """Save covariance estimates to a NetCDF file.

        Parameters
        ----------
        predictor : Predictor
            The predictor object containing binning information.
        grid : Grid
            The grid object containing binning information.
        results : CovarianceEstimates
            The covariance estimates, standard deviations, correlations, and biases.
        obs_names : list[str]
            A list of observation names.
        """
        def _sanitize(name: str) -> str:
            # Replace non-alphanumeric/underscore with underscore for NetCDF var names
            return re.sub(r"[^0-9A-Za-z_]", "_", name)

        if len(obs_names) == 1:
            obs_a, = obs_names
            obs_b = obs_a
        else:
            obs_a, obs_b = obs_names
        fname = f'{self.cov_path}_{obs_a}-{obs_b}.nc'
        dataset = netCDF4.Dataset(fname,'w') # pylint: disable=no-member
        # Global Attributes
        dataset.description = f'Dataset containing covariance estimates binned by {predictor.name}'
        dataset.history = f'Created {time.ctime(time.time())}'

        # Dimensions
        dataset.createDimension('n_grid', grid.n_gp)
        dataset.createDimension('n_predictors', predictor.n_classes)

        # Grid bins
        bins = dataset.createVariable('bins', np.float32, ('n_grid',))
        bins.setncatts({'long_name': 'Grid of correlation length', 'units': 'km'})
        bins[:] = grid.d

        # Estimates (covariance variants)
        for key, arr in results.estimates.items():
            var_name = _sanitize(f"cov_{key}")
            nc_var = dataset.createVariable(var_name, np.float32, ('n_predictors','n_grid'))
            nc_var.setncatts({'long_name': f'Covariance estimate {key}', 'units': 'none'})
            nc_var[:] = arr

        # Standard deviations (per predictor state)
        for key, arr in results.std_dev.items():
            var_name = _sanitize(f"std_dev_{key}")
            nc_var = dataset.createVariable(var_name, np.float32, ('n_predictors',))
            nc_var.setncatts({'long_name': f'Standard deviation for {key}', 'units': 'none'})
            nc_var[:] = arr

        # Correlations
        for key, arr in results.correlations.items():
            var_name = _sanitize(f"corr_{key}")
            nc_var = dataset.createVariable(var_name, np.float32, ('n_predictors','n_grid'))
            nc_var.setncatts({'long_name': f'Correlation for {key}', 'units': 'none'})
            nc_var[:] = arr

        # Bias terms (per predictor state, typically zero-separation)
        for key, arr in results.bias.items():
            var_name = _sanitize(f"bias_{key}")
            nc_var = dataset.createVariable(var_name, np.float32, ('n_predictors',))
            nc_var.setncatts({'long_name': f'Bias term for {key}', 'units': 'none'})
            nc_var[:] = arr

        dataset.close()
