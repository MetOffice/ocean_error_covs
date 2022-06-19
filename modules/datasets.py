import shutil
import os
import xarray as xr
import numpy as np
import tempfile
from timeit import default_timer as timer
from datetime import timedelta
from modules.errorCovs import ForecastErrorCovs

# Initialising the classes

class StatsDataset():
    """
    Class with methods to set up and increment an xarray dataset
    containing statistics
    """
    def __init__(self, variable, grid_shape, depthvar, num_cross_covs=0,
                 vert=False, CHUNK_STATS_IO={'y': 1000, 'x': 1000},
                 time_dim='time_counter', time_var='time_instant',
                 latname='nav_lat', lonname='nav_lon', 
                 exclude_compress = ['directions'], wrap=True,
                 PARALLEL_IO=True, COMPRESS_OUTPUT=True):
        self.variable = variable
        self.grid_shape = grid_shape
        self.num_cross_covs = num_cross_covs
        self.vert = vert
        self.depthvar = depthvar
        self.time_dim = time_dim
        self.time_var = time_var
        self.latname = latname
        self.lonname = lonname
        self.exclude_compress = exclude_compress
        self.wrap = wrap
        self.PARALLEL_IO = PARALLEL_IO
        self.CHUNK_STATS_IO = CHUNK_STATS_IO
        self.COMPRESS_OUTPUT = COMPRESS_OUTPUT
        self.list_of_vars = ["variable", latname, lonname,
                             "num_times", "times", "sum", "sumsq"]
        self.covs_ds = None
        self.covs_shape = None
        if not vert:
            prefix = "xy"
            self.dims = ["y", "x"]
            self.dircn_dims = ["N", "S", "E", "W"]
            self.covs_dims = ["dircn", "dstn", "y", "x"]
        else:
            prefix = "z"
            self.dims = [depthvar, "y", "x"]
            self.dircn_dims = ["up", "down"]
            self.covs_dims = ["dircn", "dstn", depthvar, "y", "x"]
            
        self.nam_sumsq_var = "{}_sumsq".format(prefix)
        self.nam_dists_var = "{}_dists".format(prefix)
        self.nam_covar_var = "{}_covar".format(prefix)
        self.nam_corln_var = "{}_corln".format(prefix)
        if self.num_cross_covs != 0:
            self.list_of_vars += ["directions", self.nam_sumsq_var]
            self.list_of_full_vars = self.list_of_vars + ["mean", "var", \
                                     self.nam_sumsq_var, self.nam_dists_var,
                                     self.nam_covar_var, self.nam_corln_var]


    def set_up(self, in_times):
        """
        Set up a stats xarray dataset from scratch (rather than reading it in).

        ****** PARAMETERS *****
        1. in_times:    Times of the input model forecast error fields.

        ******* RETURNS *******
        1. An xarray dataset over multiple files for the
          selected variable and level.
        """
        sum_xr = xr.DataArray(np.ma.zeros(self.grid_shape, dtype=np.float64), dims=self.dims)
        sumsq_xr = xr.DataArray(np.ma.zeros(self.grid_shape, dtype=np.float64), dims=self.dims)
        ds_dict = {"variable": self.variable,
                   "num_times": 0,
                   "times": xr.DataArray(in_times, dims=["time_counter"]),
                   "sum": sum_xr,
                   "sumsq": sumsq_xr}
        if self.num_cross_covs != 0:
            self.covs_shape = ([len(self.dircn_dims), self.num_cross_covs] + list(self.grid_shape))
            dircn_xr = xr.DataArray(self.dircn_dims, dims=["dircn"])
            sumsq_xr = xr.DataArray(np.ma.zeros(self.covs_shape, dtype=np.float64),
                                    dims=self.covs_dims)
            ds_dict.update({"directions": dircn_xr,
                            self.nam_sumsq_var: sumsq_xr})
        self.covs_ds = xr.Dataset(ds_dict).chunk(self.CHUNK_STATS_IO)


    def add_field(self, field):
        """
        Update statistics with addition of new field

        ****** PARAMETERS *****
        1. field: time slice of the input model forecast error data
        """
        self.covs_ds["num_times"] += 1
        self.covs_ds["sum"] += field
        self.covs_ds["sumsq"] += np.ma.multiply(field, field)

        if 'dstn' in self.covs_ds.dims:
            errorCovs = ForecastErrorCovs(wrap=self.wrap)
            distances = range(self.covs_ds.dims['dstn'])
            self.covs_ds[self.nam_sumsq_var].load()
            if not self.vert:
                self.covs_ds[self.nam_sumsq_var][:] += errorCovs.calc_xy_sq(self.covs_shape,
                                                                    distances, field.values)
            else:
                self.covs_ds[self.nam_sumsq_var][:] += errorCovs.calc_xyz_sq_vert(self.covs_shape,
                                                                          distances, field.values)


    def combine_stats(self, self2):
        """
        Combine statistics information from two stats datasets

        ****** PARAMETERS *****
        1. self2: dataset to be added
        """
        if self.covs_ds["variable"] != self2.covs_ds["variable"]:
            raise ValueError("Variable names in the two datasets are not the same")

        self.covs_ds["num_times"] += self2.covs_ds["num_times"]
        self.covs_ds["sum"] += self2.covs_ds["sum"]
        self.covs_ds["sumsq"] += self2.covs_ds["sumsq"]
        if 'dstn' in self.covs_ds.dims:
            if self.covs_ds.dims['dstn'] != self2.covs_ds.dims['dstn']:
                raise ValueError("Number of distances in the two datasets are not the same")
            self.covs_ds[self.nam_sumsq_var] += self2.covs_ds[self.nam_sumsq_var]


    def define_full_ds(self):
        """
        Include the full set of statistics in the xarray dataset including
        mean, variance, covariances and correlations
        """
        start_time = timer()
        errorCovs = ForecastErrorCovs(wrap=self.wrap)
        mean = errorCovs.calc_mean(self.covs_ds["num_times"], self.covs_ds["sum"])
        var = errorCovs.calc_var(self.covs_ds["num_times"], mean, self.covs_ds["sumsq"])
        end_time = timer()
        ds2_dict = {"mean": (self.dims, mean),
                    "var":  (self.dims, var)}
        print("Finished calculating the mean and variance", flush=True)
        print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)), flush=True)

        if self.num_cross_covs != 0:
            start_time = timer()
            self.covs_ds[self.nam_sumsq_var].load()
            end_time = timer()
            print("Finished loading xy_sumsq", flush=True)
            print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)), flush=True)

            start_time = timer()
            if not self.vert:
                dists = errorCovs.calc_distances(self.covs_ds[self.latname], self.covs_ds[self.lonname],
                                                 self.covs_ds[self.nam_sumsq_var])
                covs, cors = errorCovs.calc_covs_cors(self.covs_ds["num_times"], mean, var,
                                                      self.covs_ds[self.nam_sumsq_var])
            else:
                dists = errorCovs.calc_distances_vert(self.covs_ds[self.depthvar],
                                                      self.covs_ds[self.nam_sumsq_var])
                covs, cors = errorCovs.calc_covs_cors_vert(self.covs_ds["num_times"], mean, var,
                                                           self.covs_ds[self.nam_sumsq_var])
            end_time = timer()
            print("Finished calculating the covariances/correlations", flush=True)
            print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)), flush=True)

            ds2_dict.update({self.nam_dists_var: (self.covs_dims, dists),
                             self.nam_covar_var: (self.covs_dims, covs),
                             self.nam_corln_var: (self.covs_dims, cors)})

        ds2 = xr.Dataset(ds2_dict)

        return xr.merge([self.covs_ds, ds2])


    def save_stats(self, filename, save_full=False, overwrite=True):
        """
        Save data from a stats xarray dataset into a netcdf file

        ****** PARAMETERS *****
        1. filename:         Name of file to save to
        2. save_full:        Save all statistics including the mean,
                             variance and covariance
        3. overwrite:        overwrite already existing file
        """
        output_ds = self.covs_ds
        if save_full and self.nam_covar_var not in self.covs_ds.data_vars:
           output_ds = self.define_full_ds()
        self.covs_ds.close()

        if self.COMPRESS_OUTPUT:
            compress_vars = self.list_of_vars
            if save_full and self.num_cross_covs != 0:
                compress_vars = self.list_of_full_vars
            encoding = {}
            for var in compress_vars:
                if not var in self.exclude_compress:
                   encoding.update({var: {"zlib": True, "complevel": 1}})
        else:
            encoding = None

        print("Writing stats data to file {}".format(filename), flush=True)
        if overwrite:
            temp_filename = tempfile.mktemp()
            output_ds.to_netcdf(temp_filename, unlimited_dims=["time_counter"],
                                encoding=encoding)
            shutil.move(temp_filename, filename)
        else:
            output_ds.to_netcdf(filename, unlimited_dims=["time_counter"],
                                encoding=encoding)


    def read_stats(self, filename, new_times=None, read_full=False):
        """
        Read stats data from a netcdf file (produced using save_stats) into
        an xarray dataset

        ****** PARAMETERS *****
            1. filename:       Name of file to read the stats from
            2. num_cross_covs: Number of cross covariances which are going
                               to be calculated
            3. new_times:      List of times which are planned to be added to the
                               existing stats dataset
            4. read_full:      Read in the full covariance information as opposed
                               to the basic stats needed for accumulation.
        """
        input_ds = xr.open_mfdataset(filename, parallel=self.PARALLEL_IO,
                                     chunks=self.CHUNK_STATS_IO)
        vars_to_read = self.list_of_vars
        if read_full:
            vars_to_read = self.list_of_full_vars
        covs_ds = input_ds[vars_to_read]
        
        if self.num_cross_covs != 0 and covs_ds.dims['dstn'] != self.num_cross_covs:
            raise ValueError("Number of separation distances in the input ",
                             "stats file {} ".format(filename),
                             "is {} ".format(covs_ds.dims['dstn']),
                             "which is different to the number requested ",
                             "in this calculation {}".format(self.num_cross_covs))

        if self.variable not in covs_ds["variable"]:
            raise ValueError("Variable name in the ",
                             "input dataset {} ".format(self.variable),
                             "is not in the existing ",
                             "stats file {}".format(covs_ds["variable"]))

        if new_times is not None:
            if np.intersect1d(covs_ds["times"].values, new_times.values) is not None:
                raise ValueError("Some of the times in the input model error ",
                                 "fields {} ".format(new_times.values),
                                 "already exist in the output stats ",
                                 "file {}".format(covs_ds["times"].values))
            new_times_ds = xr.Dataset({"times": xr.DataArray(new_times, dims=[self.time_dim])})
            self.covs_ds = xr.concat([covs_ds, new_times_ds],
                                     dim=self.time_dim,
                                     coords=[self.time_varname],
                                     data_vars=["times"])
        else:
            self.covs_ds = covs_ds
