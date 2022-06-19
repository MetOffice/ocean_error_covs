# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
import modules.arrays as arrays
import modules.polygon as check_polygon
from modules.io_data import IO_netCDF4
from modules.utils import Utils
from time import time

# Initialising the classes
IO = IO_netCDF4()
Utils = Utils()

class CalcREs():
      """ Class of functions to calculate the representation errors """

      def __init__(self):
          self.undef = 1e10


      def calc_RE_unresolved_scales(self, args):
          """ Routine to run on multiprocessors that reads in feedback
              file and calculates representation error due to unresolved
              scales in the model

              ********* PARAMETERS ********
              1. args: dictionary with members: list of files,
                       model file, obs_type and source_type

              ********* RETURNS ***********
              1. REstats: REstats object containing the RE for each
                          model grid cell
          """
          # Read model file
          ncdata = IO.ncread_variables(args["model_file"], [args["latcor"], \
                                       args["loncor"]], dep_lev=0)
          latcor = ncdata[0]
          loncor = ncdata[1]

          if args["mask"] is None:
             mask = np.ones(latcor.shape[0], loncor.shape[1])
          else:
             mask = IO.ncread_variables(args["model_file"], [args["mask"]],
                                        dep_lev=0)[0]

          # loop over files and calculate the representation errors
          latobs = np.array([])
          lonobs = np.array([])
          obs_vals = np.array([])
          idx = args["thinning"]
          for f in args["list_of_files"]:
              print("Opening file: {}".format(f))

              # Read fdbk variables
              fdbk_var_array, depths = IO.ncread_fdbk_vars(f, args['obs_type'],
                                          args['source_types'], args['qc_val'])

              if len(fdbk_var_array.lats) > 0:

                 # Need to squash obs and model arrays to be 1D array
                 fdbk_var_array.mod_vals = fdbk_var_array.mod_vals.flatten()
                 fdbk_var_array.obs_vals = fdbk_var_array.obs_vals.flatten()
                 latobs = np.append(latobs, fdbk_var_array.lats)
                 lonobs = np.append(lonobs, fdbk_var_array.lons)
                 obs_vals = np.append(obs_vals, fdbk_var_array.obs_vals)

          if len(latobs) == 0:
             return arrays.REstats((latcor.shape[0], loncor.shape[1]))

          # Making pairs of lon, lat for obs locations
          obs_lonlat = np.array([list(pair) for pair in zip(lonobs, latobs)])

          # Subsampling obs data
          obs_lonlat = obs_lonlat[::idx,:]
          obs_vals = obs_vals[::idx]

          print("Calculating RE within a time window of {} day(s)".format(len(args["list_of_files"])))
          start_time = time()
          RE = self.calc_RE_sub_gridscale_var(latcor, loncor, mask, obs_lonlat, obs_vals,
                             args['min_num_obs'], lon_discontinuity=args['lon_discontinuity'])
          print("Calculation finished after {} mins".format((time() - start_time)/60))

          return RE


      def calc_RE_sub_gridscale_var(self, latcor, loncor, mask, obs_lonlat, obs_vals,
                                 min_num_obs, lon_discontinuity=True):
          """ Routine to calculate representation error associated with
              sub-grid scale variability

              ********* PARAMETERS ********
              1. latcor: array of latitudes at the corners of each model grid cell
              2. loncor: array of longitudes at the corners of each model grid cell
              3. mask: mask array
              4. obs_lonlat: pair of longitude and latitude for each obs
              5. obs_vals: observation values
              6. min_num_obs: min number of obs within grid cell to calculate RE
              7. lon_discontinuity: True in case there is a discontinuity in longitude

              ********* RETURNS ***********
              1. RE: RE object containing the RE and number of obs for each
                     model grid cell
          """
          idx1 = -1
          idx2 = -1

          # Initialising arrays
          RE = arrays.REstats((latcor.shape[0], loncor.shape[1]))
          for j in range(1, latcor.shape[0]):
              if lon_discontinuity:
                  idx1, idx2 = Utils.map_lon_discontinuity(loncor, j)
              for i in range(1, loncor.shape[1]):
                  if mask[j,i] == 0:
                     RE.count[j,i] = -1
                     continue

                  # Getting the grid corners
                  polygon1, polygon2 = Utils.get_grid_corners(loncor, latcor, i, j, [idx1, idx2])

                  # Determine which obs fall within the model grid cell
                  inside_cell = check_polygon.is_inside_sm_parallel(obs_lonlat, polygon1)
                  idx_obs = np.where(inside_cell)
                  obs_count = np.count_nonzero(inside_cell)
                  if polygon2 is not None:
                     inside_cell = check_polygon.is_inside_sm_parallel(obs_lonlat, polygon2)
                     idx_obs = np.append(idx_obs, np.where(inside_cell))
                     obs_count = obs_count + np.count_nonzero(inside_cell)

                  # Skipping if number of obs within a grid cell is lower than threshold
                  if obs_count < min_num_obs:
                     continue

                  # Selecting obs and estimating sub-grid scale variability
                  RE.RE[j,i] = np.nanstd(obs_vals[idx_obs])
                  RE.num_obs_in_grid[j,i] = obs_count
                  RE.count[j,i] = 1

          return RE
