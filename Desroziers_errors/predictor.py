"""Module for defining predictors used in error covariance modeling.

Author: Y Chen, University of Reading, 2025
"""
import configparser
from collections import OrderedDict
import typing

import numpy as np

from . import log

class Predictor:
    """Background and observation error covariance will depend on state of the
    system. Different regimes of weather patterns or different ocean conditions
    may have different error characteristics. For example, the errors can depend
    on ocean bathemetry or optical water type (OWT).

    For continuous predictors (e.g. OWT), their values will be binned into
    discrete ranges. For categorical predictors, the categories will be used
    directly.
    """
    def __init__(self, configs: dict[str, configparser.ConfigParser], mask_func=None) -> None:
        self.max_obs_per_day = OrderedDict()
        n_classes: list[int] = []
        predictor_names: list[str] = []
        for name, config in configs.items():
            predictor_cfg = config['Predictor']
            predictor_names.append( predictor_cfg.get('name', 'predictor') )
            n_classes.append(predictor_cfg.getint('n_classes', 5))
            self.max_obs_per_day[name] = predictor_cfg.getint('max_obs_per_day', 2000)

        assert len(set(n_classes)) == 1, \
                "All observations must have the same number of classes for predictors."

        assert len(set(predictor_names)) == 1, \
            "All observations must have the same predictor name."

        self.name = predictor_names[0]
        self.n_classes = n_classes[0]
        self._rng = np.random.default_rng()
        self.mask_func = mask_func if mask_func is not None else self.get_mask

    def get_mask(self, _name:str, data:dict[str, np.ndarray], i:int) -> np.ndarray:
        """Default mask function for i-th predictor bin.

        This function defines a mask that selects all observations.

        Parameters
        ----------
        name: str
            Name of the observation type. This is specified in the
            .ini file for each observation type.
        _data : dict[str, np.ndarray]
            Dictionary containing all data variables from input handler.
        grid_t : Grid
            Grid object containing observation metadata.
        _i : int
            i-th predictor bin.

        Returns
        -------
        np.ndarray
            Boolean array for selecting observations in the i-th predictor bin.
        """
        if 'predictor'in data:
            return data['predictor'] == i
        else:
            return np.ones(len(data['d_ob_j']), dtype=bool)

    def get_sample_indices(
        self, i_day: int, data: dict[str, dict[str, np.ndarray]]
    ) -> typing.Iterator[dict[str, np.ndarray]]:
        """Iterator that yields indices of samples for each predictor category.

        This method iterates through each class of the predictor and returns the indices
        of observations that fall within each class. If the number of observations in a
        class exceeds the maximum allowed per day, the observations are randomly shuffled
        and limited to the maximum threshold.

        Parameters
        ----------
        i_day : int
            Day index for logging purposes.
        data : dict[str, np.ndarray]
            Collection containing at least the ``predictor`` array for each
            observation type.

        Yields
        ------
        dict[str, np.ndarray]
            Each entry of dictionary contains an array of indices for given
            observation type belonging to the current predictor class.
            The array may be limited to max_obs_per_day elements
            if the class contains too many observations.
        """
        for idx in range(self.n_classes):
            log.logger.info(f'Starting sampling for the {idx + 1}-th predictor class')

            sampled_indices: dict[str, np.ndarray] = OrderedDict()

            for name, max_obs in self.max_obs_per_day.items():
                # get mask for current predictor class
                mask = self.mask_func(name, data[name], idx)
                # get sampled indices
                sampled_indices[name] = np.nonzero(mask)[0]
                n_obs_this_day = len(sampled_indices[name])
                # limiting the number of observations per file
                if n_obs_this_day <= max_obs:
                    log.logger.info(
                        f'Number of {name} obs in {idx + 1}-th predictor class '
                        f'for {i_day}-th file: {n_obs_this_day}'
                    )
                    continue
                log.logger.warning(
                    f'Too many {name} obs. in {idx + 1}-th predictor class for '
                    f'{i_day}-th file: {n_obs_this_day}, limiting to {max_obs}'
                )
                sampled_indices[name] = self._rng.choice(sampled_indices[name],
                                                         size=max_obs,
                                                         replace=False
                                                         )
                log.logger.info(
                    f'Number of {name} obs in {idx + 1}-th predictor class '
                    f'for {i_day}-th file: {max_obs}'
                )

            yield sampled_indices
