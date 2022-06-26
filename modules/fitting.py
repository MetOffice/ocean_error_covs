# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
import numpy.ma as ma
import scipy.optimize as opt
import xarray as xr
import modules.functions as functions
from modules.io_data import IO_xarray

class FittingHL():
      """ Class of fitting functions for HL method"""

      def __init__(self):
         self.MASK_VAL = -1e10  


      def calc_x_positions(self, bins):
          """ Calculate x positions based on the bins of separation distance

          ******** PARAMETERS **********
          1. bins: bins of separation distance

          ******** RETURNS *************
          1. x_val: values of x positions
          """
          x_val=[0.5*bins[0]]
          for n in range(0, len(bins)-1):
              x_val += [0.5 * (bins[n] + bins[n+1])]
          x_val = np.array(x_val)
          return x_val


      def create_arg_list(self, x_val, cors, var, numobsvar, min_num_obs, 
                          func_name, num_funcs, lenscale, max_iter, scalefac):
          """ Create list of arguments for running in parallel

          ******** PARAMETERS **********
          1. xval: list of x positions
          2. cors: correlation of innovation statistics
          3. var: grid box variance
          4. numobsvar: number of obs within grid box
          5. min_num_obs: minimum number of obs to carry on with fitting
          6. func_name: name of fitting function
          7. num_funcs: number of functions to use
          8. lenscale: lengthscales to apply
          9. max_iter: max number of iterations
          10. scalefac: factor to scale the variances in case it is needed

          ******** RETURNS *************
          1. arg_list: list with arguments
          """
          arg_lists = [] 
          for nn in range(0, cors.shape[0]):
              for pp in range(0, cors.shape[1]):

                  if numobsvar[nn,pp] < min_num_obs:
                      cors[nn,pp,:] = np.nan
                      var[nn,pp] = np.nan

                  arg_lists += [{ "covs": cors[nn,pp,:]*var[nn,pp]*scalefac, # NOTE: Reweighting by central point
                                  "x_vals": x_val[:],
                                  "var": var[nn,pp]*scalefac,
                                  "max_iter": max_iter }]

                  lscales = []
                  for nfunc in range(0, num_funcs):
                      lscales += [lenscale[nfunc][nn,pp]] if isinstance(lenscale[nfunc], np.ndarray) \
                                                         else [lenscale[nfunc]]

                  if func_name == 'MultiGauss':
                      i_guess = []
                      for i in range(0, num_funcs):
                          i_guess.append(0.5*var[nn,pp]*scalefac)
                          i_guess.append(lscales[i])
                      arg_lists[-1]["func"] = functions.MultiGaussFunction(initial_guess=i_guess, 
                                                              weights=None, num_funcs=num_funcs)
                  elif func_name == 'MultiGauss_Fixed':
                      i_guess = [0.5*var[nn,pp]*scalefac] * num_funcs
                      arg_lists[-1]["func"] = functions.MultiGaussFunction_FixedLenScale(initial_guess=i_guess, 
                                                        weights=None, num_funcs=num_funcs, lenscales=lscales)
                  else:
                    raise ValueError('[ERROR] FUNCTION NAME SELECTED IS NOT AVAILABLE')

          return arg_lists


      def fitter(self, args):
          """
          Algorithm to do the fitting

          ************ PARAMETERS ***********
          1. args: dictionary containing correlations and variance for 1 point,
                   instance of fitting function and x values

          ************* RETURNS *************
          1. fit: fitting results
          2. obs_err: obs error 
          3. rss_vs_func: Residual sum of squares from fitting
          4. rss_vs_mean: Residual sum of squares from mean covs
          """
          
          # Read arguments
          covs = args["covs"]
          func = args["func"]
          x_vals = args["x_vals"]
          var = args["var"]
          max_iter = args["max_iter"]
              
          #check the corrs for nans and calculate fit
          nans = np.isnan(covs)
          num_points = np.sum(np.logical_not(nans))
          if num_points < func.num_params() + 1:
              fit = self.MASK_VAL*np.ones((func.num_params()))
              obs_err = self.MASK_VAL
              rss_vs_func  = self.MASK_VAL
              rss_vs_mean = self.MASK_VAL
              num_points = int(self.MASK_VAL)
          else:
              #NOTE: we do not use the zero point in the fitting in order to get the obs error
              fit = opt.fmin_tnc(func.cost_func, func.init_guess(), disp=0, fprime=func.cost_grad,
                                 maxfun=max_iter, args=(x_vals[:], covs[:]),
                                 bounds=func.bounds())[0]
              obs_err = var - func.func(0, *fit)

              # Calculate Residual Sum of Squares (RSS) values.
              # Note: for the RSS of the function this is the same as the
              # chi_sq statistic if all 'variances' are 1.

              rss_vs_func = func.chi_statistic(fit, x_vals[:], 1., covs[:])
              rss_vs_mean = np.sum((covs - np.mean(covs))**2.)

          return (fit, obs_err, rss_vs_func, rss_vs_mean, num_points)


      def results_to_grid(self, results, nlat, nlon):
          """ Grid the results from multiprocessing 

          ************ PARAMETERS ***********
          1. results: results merged over all the processors
          2. nlat: number of latitudes
          3. nlon: number of longitudes

          ************* RETURNS *************
          1. params: fitting parameters gridded
          2. obs_err: obs error gridded
          3. chi: chi squares gridded
          4. dof: degrees of freedom for the fit
                  ( = Number of points - number of parameters)
          """

          rss_func_grid = ma.zeros((nlat, nlon))
          rss_mean_grid = ma.zeros((nlat, nlon))
          obs_err=ma.zeros((nlat, nlon))
          dof_grid = ma.zeros((nlat, nlon), dtype=int)
          num_params = len(results[0][0][:])
          params = []
          for param in range(0, num_params):
              params += [ma.zeros((nlat, nlon))]

          n = 0
          for nn in range(0, nlat):
              for pp in range(0, nlon):
                for param in range(0, num_params):
                    params[param][nn,pp] = results[n][0][param]
                obs_err[nn,pp] = results[n][1]
                rss_func_grid[nn,pp] = results[n][2]
                rss_mean_grid[nn,pp] = results[n][3]
                dof_grid[nn,pp] = results[n][4] - num_params
                n += 1

          return params, obs_err, rss_func_grid, rss_mean_grid, dof_grid


class FittingForecastDiff():
      """ Class of fitting functions for methods based on forecast differences"""

      def __init__(self, method="fmin_tnc", func_name='MultiGauss',
                   anisotropic=True, num_funcs=2, lenscale=(40, 400),
                   vert=False, max_iter=1000, min_valid_points=10,
                   level=None, scalefac=1.0):
          self.method = method
          self.func_name = func_name
          self.anisotropic = anisotropic
          self.num_funcs = num_funcs
          self.lenscale = lenscale
          self.vert = vert
          self.max_iter = max_iter
          self.min_valid_points = min_valid_points
          self.level = level
          self.scalefac = scalefac
          if self.anisotropic and self.vert:
             self.directions = ["Up", "Down"]
             self.varout = ["params_up", "params_down"]
          elif self.anisotropic:
             self.directions = ["W-E", "S-N"]
             self.varout = ["params_x", "params_y"]
          else:
             self.directions = ["all"]
             self.varout = ["params_iso"]
          varnames = []
          for v in range(0, num_funcs):
              if func_name == "MultiGauss_Fixed":
                 varnames.extend(["var{}".format(v+1)])
              else:
                 varnames.extend(["var{}".format(v+1), "len{}".format(v+1)])
          self.varnames = varnames


      def select_vars(self, covs_ds):
          """ Selecting variables from dataset for function fitting.

          ******** PARAMETERS **********
          1. covs_ds: xarray dataset

          ******** RETURNS *************
          1. covs: covariances
          2. xvals: distances
          3. lats: latitudes
          4. lons: longitudes
          """
          if self.vert:
              if self.level is None:
                  raise ValueError("[ERROR] For vertical covariance a model level is needded")
              io = IO_xarray()
              varlist = list(covs_ds.variables.keys())
              z = io.findvardepth(varlist)
              indexers = {z: self.level}
              covs = covs_ds["z_covar"].isel(**indexers)
              xvals = covs_ds["z_dists"].isel(**indexers)
          else:
              covs = covs_ds["xy_covar"]
              xvals = covs_ds["xy_dists"]

          return np.array(covs), np.array(xvals), covs_ds["nav_lat"], covs_ds["nav_lon"]


      def average_covs(self, xvals, covs):
          """ Averaging covariances for each or all directions.

          ******** PARAMETERS **********
          1. xvals: x positions
          2. covs: covariances

          ******** RETURNS *************
          1. xvals: distances averaged for each or all directions
          2. covs: covariances averaged for each or all directions
          """
          covs_avg = []
          xvals_avg = []
          if self.anisotropic and self.vert:
             xvals_avg.append(xvals[0,:])
             covs_avg.append(covs[0,:]*self.scalefac)
             xvals_avg.append(xvals[1,:])
             covs_avg.append(covs[1,:]*self.scalefac)
          elif self.anisotropic:
             xvals_avg.append(np.ma.average(xvals[2:4,:], axis=0))
             covs_avg.append(np.ma.average(covs[2:4,:], axis=0)*self.scalefac)
             xvals_avg.append(np.ma.average(xvals[0:2,:], axis=0))
             covs_avg.append(np.ma.average(covs[0:2,:], axis=0)*self.scalefac)
          else:
             xvals_avg.append(np.ma.average(xvals, axis=0))
             covs_avg.append(np.ma.average(covs, axis=0)*self.scalefac)

          return xvals_avg, covs_avg


      def create_arg_list(self, xvals, covs):
          """ Create list of arguments for running function
              fitting in parallel.

          ******** PARAMETERS **********
          1. xvals: distances
          2. covs: covariances

          ******** RETURNS *************
          1. arg_list: argument list with all parameters
                       for function fitting
          """
          arg_lists = []
          for nn in range(0, covs.shape[1]):
              for pp in range(0, covs.shape[2]):
                  arg_lists += [{ "covs": covs[:,nn,pp],
                                  "x_vals": xvals[:,nn,pp],
                                  "max_iter": self.max_iter,
                                  "method": self.method,
                                  "min_valid_points": self.min_valid_points}]

                  lscales = []
                  for nfunc in range(0, self.num_funcs):
                      lscales += [self.lenscale[nfunc][nn,pp]] if \
                                 isinstance(self.lenscale[nfunc], np.ndarray) else \
                                 [self.lenscale[nfunc]]

                  if self.func_name == 'MultiGauss':
                      i_guess = []
                      for i in range(0, self.num_funcs):
                          i_guess.extend([0.5*self.scalefac, lscales[i]])
                      arg_lists[-1]["func"] = functions.MultiGaussFunction(initial_guess=i_guess,
                                                          weights=None, num_funcs=self.num_funcs)
                  elif self.func_name == 'MultiGauss_Fixed':
                      i_guess = [0.5*self.scalefac]*self.num_funcs
                      arg_lists[-1]["func"] = functions.MultiGaussFunction_FixedLenScale(initial_guess=i_guess,
                                                     weights=None, num_funcs=self.num_funcs, lenscales=lscales)
                  elif self.func_name == 'Gauss':
                       i_guess = [0.5*self.scalefac, lscales[0]]
                       arg_lists[-1]["func"] = functions.GaussFunction(initial_guess=i_guess, weights=None)
                  else:
                    raise ValueError('[ERROR] FUNCTION NAME SELECTED IS NOT AVAILABLE')

          return arg_lists


      @staticmethod
      def fitter(args):
          """
          Algorithm to do the fitting

          ************ PARAMETERS ***********
          1. args: dictionary containing stats for 1 point,
                   instance of fitting function and x values

          ************* RETURNS *************
          1. fit: fitting results
          2. chi: chi values
          """
          func = args["func"]
          min_valid_points = args["min_valid_points"]
          max_iter = args["max_iter"]
          method = args["method"]
          covs_v = np.array(args["covs"])
          xvals_v = np.array(args["x_vals"])

          # Check the corrs for nans and calculate fit
          nans = np.isnan(covs_v)
          if (len(covs_v) - np.sum(nans)) < min_valid_points:
              fit = np.nan*np.ones((func.num_params()))
              chi = np.nan
          else:
              xvals_v = xvals_v[~nans]
              covs_v = covs_v[~nans]

              if method == "fmin_tnc":
                 fit = opt.fmin_tnc(func.cost_func, func.init_guess(), disp=0, fprime=func.cost_grad,
                                    maxfun=max_iter, args=(xvals_v, covs_v), bounds=func.bounds())[0]
              elif method == "L-BFGS-B":
                 fit = opt.minimize(func.cost_func, func.init_guess(), method=method, jac=func.cost_grad,
                                    args=(xvals_v, covs_v), bounds=func.bounds())['x']
              else:
                 raise ValueError('[ERROR] MINIMISER METHOD NOT FOUND!')

              chi = func.cost_func(fit, xvals_v[:], covs_v[:])

          return (fit, chi)


      @staticmethod
      def results_to_grid(results, nlat, nlon):
          """ Grid the results from multiprocessing

          ************ PARAMETERS ***********
          1. results: results merged over all the processors
          2. nlat: number of latitudes
          3. nlon: number of longitudes

          ************* RETURNS *************
          1. params: fitting parameters gridded
          2. chi: chi values gridded
          """
          chi_grid = ma.zeros((nlat, nlon))
          num_params = len(results[0][0][:])
          params = []
          for param in range(0, num_params):
              params += [ma.zeros((nlat, nlon))]
          n = 0
          for nn in range(0, nlat):
               for pp in range(0, nlon):
                   for param in range(0, num_params):
                       params[param][nn,pp] = results[n][0][param]
                   chi_grid[nn,pp] = results[n][1]
                   n += 1

          return params, chi_grid


      @staticmethod
      def create_fitting_ds(params, lats, lons, dsname, varnames, to_merge=None):
          """ Create xarray dataset with fitting results

          ************ PARAMETERS ***********
          1. params: parameter values from function fitting
          2. lats: latitudes
          3. lons: longitudes
          4. dsname: name of the xarray dataset
          5. varnames: names of the parameters in the dataset
          6. to_merge: dataset to be merged

          ************* RETURNS *************
          1. params_out: output dataset
          """
          params_out = xr.DataArray(params, name=dsname, coords={"params": varnames, \
                         "nav_lat": lats, "nav_lon": lons}, dims=["params", "y", "x"])
          if to_merge is not None:
             params_out = xr.merge([to_merge, params_out])

          return params_out
