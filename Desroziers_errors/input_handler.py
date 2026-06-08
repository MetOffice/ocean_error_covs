"""Dealing with input and output of NetCDF files

Author: Y Chen, University of Reading, 2025
"""
import configparser
import datetime
import os
import typing

import netCDF4
import numpy as np

from . import log


class _DerivationSource(typing.TypedDict):
    d: tuple[str, str]
    sign: typing.Callable[[np.ndarray, np.ndarray], np.ndarray]


class InputHandler:
    """Class to handle input NetCDF files based on configuration.

    This class read configuration from a single configuration section,
    which corresponds to a single input data source, and fixed variables.
    """
    def __init__(self, config:configparser.SectionProxy, estimates: list[str]) -> None:
        """Initialize InputHandler with configuration.

        Parameters
        ----------
        config : configparser.SectionProxy
            Configuration section for the input data source.
        estimates : list[str]
            List of estimate types to process.
        """
        _innov = {'HBH': 'd_ab', 'R': 'd_oa'}
        d_names = [_innov[estimate] for estimate in estimates if estimate in _innov]
        d_names.append('d_ob_j')
        self.input_dir = config.get('input_dir', './input_data/')
        self.filename_format = config.get('filename_format',
                        'data_%Y%m%dT%H%MZ_Bathy_data.nc')
        self.data_freq = config.get('data_freq', '1D')  # in hours
        self.data_start = config.get('data_start', '2018-04-01T00:00:00.00')
        self.data_end = config.get('data_end', '2018-04-01T00:00:00.00')
        self.variables = config.get('variables', '').split(',')
        self.lon_name = config.get('lon_name', 'lon')
        self.lat_name = config.get('lat_name', 'lat')
        self.vert_name = config.get('vert_name', 'vert')
        self.pred_name = config.get('pred_name', None)

        # deal with the existence of only two of d_ob, d_oa, d_ab
        self.names = {}
        for vname in ['d_ob', 'd_oa', 'd_ab']:
            nameval = config.get(f'{vname}_name', 'none')
            if nameval != 'none':
                if vname != 'd_ob':
                    self.names[vname] = nameval
                else:
                    self.names[vname+'_j'] = nameval

        derivation_sources: dict[str, _DerivationSource] = {
            'd_oa': {'d': ('d_ob_j', 'd_ab'), 'sign': np.subtract},
            'd_ab': {'d': ('d_ob_j', 'd_oa'), 'sign': np.subtract},
            'd_ob_j': {'d': ('d_oa', 'd_ab'), 'sign': np.add}
        }

        self._innov_to_read = [vname for vname in d_names if vname in self.names]
        self._innov_to_derive: dict[str, _DerivationSource] = {vname: derivation_sources[vname]
                                 for vname in d_names
                                 if vname not in self.names}
        self.other_vars = [var for var in self.variables
                           if var not in self._innov_to_read and
                              var not in [self.lon_name, self.lat_name,
                                          self.vert_name, self.pred_name]
                           ]
        # Precompute innovation handling once to keep per-file read loop simple.
        self.fnames = self.get_fnames()


    def get_fnames(self) -> typing.Iterator[str]:
        """Generate list of filenames based on configuration.

        Yields
        ------
        str
            List of generated filenames.
        """
        if self.data_freq == '':
            for f in os.listdir(self.input_dir):
                yield os.path.join(self.input_dir, f)
        else:
            freq = typing.cast(typing.Literal['h', 'D', 'M', 'Y'], self.data_freq[-1])
            start_dt = np.datetime64(self.data_start).astype(f'datetime64[{freq}]')
            end_dt = np.datetime64(self.data_end).astype(f'datetime64[{freq}]')
            delta = np.timedelta64(int(self.data_freq[:-1]), freq)
            all_dates = typing.cast(typing.Iterator[np.datetime64],
                                      np.arange(start_dt, end_dt + delta, delta)
                                    )
            for dt in all_dates:
                yield os.path.join(
                    self.input_dir,
                    dt.astype(datetime.datetime).strftime(self.filename_format)
                )

    def read(self, is_horizontal: bool) -> typing.Iterator[dict[str, np.ndarray]]:
        """Read specified variables from a NetCDF file.

        The filenames are derived from :py:meth:`input.InputHandler.get_fnames`
        and variables are given in the config file.

        Yields
        ------
        dict[str, np.ndarray]
            Dictionary of variable names and their corresponding data arrays over
            all files.
        """
        # todo: time dimension handling?
        self.fnames = self.get_fnames()
        for fname in self.fnames:
            data = {}
            msg = f'Reading {fname}'
            log.logger.info(msg)
            with netCDF4.Dataset(fname, 'r') as dataset: # pylint: disable=no-member
                if is_horizontal:
                    data['lon'] = dataset.variables[self.lon_name][:]
                    data['lat'] = dataset.variables[self.lat_name][:]
                else:
                    data['vert'] = dataset.variables[self.vert_name][:]

                if self.pred_name is not None:
                    data['predictor'] = dataset.variables[self.pred_name][:]
                for var in self.other_vars:
                    data[var] = dataset.variables[var][:]

                for vname in self._innov_to_read:
                    item = self.names[vname]
                    data[vname] = dataset.variables[item][:]

                for vname in self._innov_to_derive:
                    src_a, src_b = self._innov_to_derive[vname]['d']
                    src_a = self.names[src_a]
                    src_b = self.names[src_b]
                    sign = self._innov_to_derive[vname]['sign']
                    data[vname] = sign( dataset.variables[src_a][:], dataset.variables[src_b][:])

            yield data


if __name__ == "__main__":
    pass  # for testing purposes only
