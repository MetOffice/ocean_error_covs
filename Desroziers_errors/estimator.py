"""
Estimator for calculating Desroziers error covariances.

Author: Y Chen, University of Reading, 2025
"""
from collections import OrderedDict
from dataclasses import dataclass
from typing import Iterator

import numpy as np

from .grid import Grid
from .predictor import Predictor


@dataclass
class BinnedStats:
    """Container for binned statistics."""
    n_sample_obs: dict[str, np.ndarray]
    n_obs_bin: np.ndarray
    err_sum: dict[str, np.ndarray]
    d_sum: dict[str, np.ndarray]

@dataclass
class CovarianceEstimates:
    """Container for covariance estimates."""
    estimates:dict[str, np.ndarray]
    std_dev: dict[str, np.ndarray]
    correlations: dict[str, np.ndarray]
    bias: dict[str, np.ndarray]


class CovStatBinner:
    """Collect statistics for desroziers diagnostics.

    This class bins raw data into predictor bins and spatial correlation grid.
    """
    def __init__(self, grid: Grid, predictor: Predictor, obs_names: list[str],
                 estimates: list[str]) -> None:
        """Initialize with only the components needed for computation.

        Parameters
        ----------
        grid : Grid
            Spatial grid for distance binning
        predictor : Predictor
            Predictor for state binning
        obs_names : list[str]
            List of observation types to process
        estimates : list[str]
            List of estimate types to process
        """
        self.grid = grid
        self.predictor = predictor
        self.obs_names = obs_names
        self.estimates = estimates
        _innov = {'HBH': 'd_ab', 'R': 'd_oa', 'R+HBH': 'd_ob_j'}
        self.innovs = [_innov[v] for v in self.estimates]


    def do_binning(self, data_iterator: Iterator[dict]) -> BinnedStats:
        """Collect error statistics required by the Desrozier diagnostics

        Parameters
        ----------
        input : string
            Path to innovation files.
        predictor : Predictor
            Predictor object used to bin data.

        Returns
        -------
        None.
        """
        # Initialize accumulator arrays
        n_state = self.predictor.n_classes
        n_grid = self.grid.n_gp

        n_sample_obs =  {name: np.zeros(n_state) for name in self.obs_names}
        n_obs_bin = np.zeros((n_state, n_grid))
        err_sum = self._init_error_dict(n_state, n_grid)
        d_sum = self._init_difference_dict(n_state, n_grid)

        # Process each day
        for day_idx, data in enumerate(data_iterator):
            self.grid.update_coord(data)

            daily_results = self._process_single_file(day_idx, data, n_sample_obs)
            # Accumulate results
            n_obs_bin += daily_results['n_obs_bin']
            self._accumulate_dicts(err_sum, daily_results['err_sum'])
            self._accumulate_dicts(d_sum, daily_results['d_sum'])

        return BinnedStats(n_sample_obs, n_obs_bin, err_sum, d_sum)

    def _process_single_file(self, day_idx: int, data: dict[str, dict[str, np.ndarray]],
                             n_sample_obs: dict[str, np.ndarray]) -> dict:
        """Process data within a single file.

        Depending on the file frequency. This could be daily/monthly, etc.
        This function bins data into predictor bins.

        Parameters
        ----------
        day_idx : int
            Day index in the data iterator
        data : dict
            Observation data for the day
        n_sample_obs : np.ndarray
            Accumulated sample observation counts (modified in place)

        Returns
        -------
        dict
            Single file results including n_obs_bin, err_sum, and d_sum
        """
        n_state = self.predictor.n_classes
        n_grid = self.grid.n_gp

        n_obs_bin_daily = np.zeros((n_state, n_grid))
        err_sum_daily = self._init_error_dict(n_state, n_grid)
        d_sum_daily = self._init_difference_dict(n_state, n_grid)

        # Process each predictor bin
        sample_iterator = self.predictor.get_sample_indices(day_idx, data)

        for state_idx, sample_indices in enumerate(sample_iterator):
            n_obs = [len(sample_indices[name]) == 0 for name in sample_indices]
            if any(n_obs):
                continue

            for name in sample_indices:
                n_sample_obs[name][state_idx] += len(sample_indices[name])

            # Compute Desroziers diagnostics
            n_obs_bin_daily[state_idx], d_sum_state, err_uc_state = \
                self._compute_desroziers(self._extract_differences(data,
                                                                   sample_indices
                                                                   ),
                                         sample_indices
                                         )
            # Store daily results
            for key in err_sum_daily:
                err_sum_daily[key][state_idx] = err_uc_state[key]
            for key in d_sum_daily:
                d_sum_daily[key][state_idx] = d_sum_state[key]

        return {'n_obs_bin': n_obs_bin_daily,
                'err_sum': err_sum_daily,
                'd_sum': d_sum_daily
        }

    def _compute_desroziers(self, d: dict, sample_indices: dict[str, np.ndarray],
                            chunk_size: int = 500) -> tuple:
        """Data binning for correlation grid, and compute un-centred expectations.

        Parameters
        ----------
        d : dict[str, np.ndarray]
            - vector of innovations, y - Hx^b, key: 'd_ob'
            - analysis innovations, y - Hx^a, key: 'd_oa'
            - increments, Hx^a - Hx^b, key: 'd_ab'
        sample_indices : np.ndarray
            Indices for observations belonging to the current predictor bin
        chunk_size : int
            Number of samples to process in each chunk to manage memory usage

        Returns
        -------
        n_sample_grid : np.ndarray
            number of samples (pairs of observations) within each distance category.
            shape: len(grid)
        d_sum : dict[str, np.ndarray]
            sums of d_ob, d_oa, d_ab within each distance bin, keys: 'd_ob', 'd_oa', 'd_ab'
            Each distance category has two points on space.
            shape: len(grid) x 2
        err_uc : dict[str, np.ndarray]
            products by un-centred samples for R, R+HBH, HBH, keys: 'R', 'R+HBH', 'HBH'
            shape: len(grid)
        """
        n_grid = self.grid.n_gp
        n_samples = [len(sample_indices[name]) for name in sample_indices]
        if len(n_samples) == 1:
            name_a, = list(sample_indices.keys())
            name_b = name_a
            n_samples.append(n_samples[0])
        else:
            (name_a, name_b) = sample_indices.keys()

        n_a, n_b = n_samples[0], n_samples[1]

        n_sample_grid = np.zeros(n_grid, dtype=np.intp)
        d_sum = {v: np.zeros(n_grid) for v in ['d_ob_i'] + self.innovs}
        err_uc = {v: np.zeros(n_grid) for v in self.estimates}

        # Pre-extract 1-D b-side vectors to avoid repeated dict lookups
        d_b = {v: d[name_b][v] for v in self.innovs}
        d_ob_i_a = d[name_a]['d_ob_j']  # shape (n_a,)

        for i0 in range(0, n_a, chunk_size):
            i1 = min(i0 + chunk_size, n_a)
            ci = i1 - i0
            chunk_indices_a = {name_a: sample_indices[name_a][i0:i1]}
            chunk_indices_b = {name_b: sample_indices[name_b]}
            chunk_bin_idx = self.grid.calc_distance(chunk_indices_a, chunk_indices_b)

            n_sample_grid += np.bincount(chunk_bin_idx, minlength=n_grid)

            # ravel() on a broadcast (non-contiguous) array returns a copy — ci*n_b floats
            d_ob_i_chunk = np.broadcast_to(d_ob_i_a[i0:i1, None], (ci, n_b)).ravel()
            np.add.at(d_sum['d_ob_i'], chunk_bin_idx, d_ob_i_chunk)

            # Compute and cache b-side chunks (reused for both d_sum and err_uc)
            d_b_chunks = {v: np.broadcast_to(d_b[v][None, :], (ci, n_b)).ravel()
                        for v in self.innovs}
            for v in self.innovs:
                np.add.at(d_sum[v], chunk_bin_idx, d_b_chunks[v])

            for v_err, d_name in zip(self.estimates, self.innovs):
                np.add.at(err_uc[v_err], chunk_bin_idx,
                        d_ob_i_chunk * d_b_chunks[d_name])

        # remove means and compute covariances
        # this can be fully vectorised
        # c = n_sample_grid
        # for v, v_err in zip(['d_oa', 'd_ab', 'd_ob_j'], ['R', 'HBH', 'R+HBH']):
        #     err[v_err] = (err_uc[v_err] - (d_sum['d_ob_i']*d_sum[v])/c)/(c - 1)

        return n_sample_grid, d_sum, err_uc

    # Helper methods
    def _init_error_dict(self, n_state: int, n_grid: int) -> dict[str, np.ndarray]:
        return {v: np.zeros((n_state, n_grid)) for v in self.estimates}

    def _init_difference_dict(self, n_state: int, n_grid: int) -> dict[str, np.ndarray]:
        return {v: np.zeros((n_state, n_grid))
                for v in ['d_ob_i'] + self.innovs
                }

    def _extract_differences(self, data: dict, indices: dict[str, np.ndarray]) -> dict:
        return { name:
                    {
                      v: data[name][v][indices[name]] for v in set(self.innovs + ['d_ob_j',])
                    }
                for name in data
                }

    @staticmethod
    def _accumulate_dicts(target: dict, source: dict) -> None:
        for key in target:
            target[key] += source[key]


class CorrEstimator:
    """Converts binned data into covariance and correlation estimates.
    """
    def __init__(self, estimates: list, ss_lim: int, smooth_window: int) -> None:
        self.estimates = estimates
        self.ss_lim = ss_lim
        self.smooth_window = smooth_window
        _innov = {'HBH': 'd_ab', 'R': 'd_oa', 'R+HBH': 'd_ob_j'}
        self.innovs = [_innov[v] for v in self.estimates]

    def do_compute(self, results: BinnedStats) -> CovarianceEstimates:
        """Compute covariance estimates from accumulated statistics.

        Removes bias and computes various error covariance estimates
        including options with/without bias correction.

        Parameters
        ----------
        results : CorrelationResults
            Accumulated correlation results

        Returns
        -------
        CovarianceEstimates
            Estimates, standard deviations, correlations, and biases
        """
        n_state, n_grid = results.n_obs_bin.shape
        est = self._init_estimate_dict(n_state, n_grid)
        sdev = {v: np.zeros((n_state)) for v in est}
        corr = {v: np.zeros((n_state, n_grid)) for v in est}
        bias = {v: np.zeros(n_state) for v in self.innovs}

        mask = results.n_obs_bin > self.ss_lim
        c = results.n_obs_bin[mask]

        # Compute bias-corrected estimates
        for v, v_err in zip(self.innovs, self.estimates):
            est[v_err][mask] = (
                results.err_sum[v_err][mask] -
                (results.d_sum['d_ob_i'][mask] * results.d_sum[v][mask]) / c
            ) / (c - 1)
            est[v_err][~mask] = np.nan
            corr[v_err][~mask] = np.nan
            sdev[v_err], corr[v_err] = self._convert_cov_to_corr(est[v_err])

        # Compute innovation bias terms
        mean_innov = np.mean(
            np.where(mask, results.d_sum['d_ob_i'] / results.n_obs_bin, 0.), axis=1
        )
        innov_bias_cov = np.zeros((n_state, n_grid))
        innov_bias_cov[:] = mean_innov[:, None] ** 2

        c = results.n_obs_bin
        # Compute variants with/without bias
        for v in set(self.estimates) & {'HBH', 'R+HBH'}:
            # No bias removed
            est[f'{v}_no_bias_removed'][mask] = results.err_sum[v][mask] / (c[mask] - 1)
            # With innovation bias
            est[f'{v}_with_innov_bias'][mask] = est[v][mask] + c[mask] * innov_bias_cov[mask] / (c[mask] - 1)

            # Convert to correlations
            for suffix in ['_with_innov_bias', '_no_bias_removed']:
                key = f'{v}{suffix}'
                est[key][~mask] = np.nan
                corr[key][~mask] = np.nan
                sdev[key], corr[key] = self._convert_cov_to_corr(est[key])

        # Compute bias at zero separation
        mask_zero = results.n_obs_bin[:, 0] > self.ss_lim
        for v in self.innovs:
            bias[v][mask_zero] = (
                results.d_sum[v][mask_zero, 0] / results.n_obs_bin[mask_zero, 0]
            )
            bias[v][~mask_zero] = np.nan

        return CovarianceEstimates(est, sdev, corr, bias)

    def _convert_cov_to_corr(self, cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Convert covariance to correlation and standard deviation.

        Parameters
        ----------
        cov : np.ndarray
            Covariance matrix

        Returns
        -------
        tuple
            (standard deviations, smoothed correlations)
        """
        var = cov[:, 0]
        sdev = np.sqrt(var)
        corr = np.minimum(cov / var[:, None], 1.)

        # Smooth correlations
        corr_sm = np.empty(cov.shape)
        corr_sm[:, 0] = 1.
        sw = self.smooth_window # number of points to smooth over, better if even
        n_state, n_grid = cov.shape
        for i in range(n_state):
            for j in range(1, n_grid):
                j0 = max(0, j - sw//2)
                j1 = min(n_grid, j + sw//2)
                corr_sm[i,j]=np.mean(corr[i, np.arange(j0, j1)])

        return sdev, corr_sm

    def _init_estimate_dict(self, n_state: int, n_grid: int) -> dict[str, np.ndarray]:
        """Initialize estimate dictionary with all variants."""
        keys = self.estimates + \
               [f'{v}_with_innov_bias'
                for v in set(self.estimates) & {'HBH', 'R+HBH'}]+ \
               [f'{v}_no_bias_removed'
                for v in set(self.estimates) & {'HBH', 'R+HBH'}]
        return {v: np.zeros((n_state, n_grid)) for v in keys}
