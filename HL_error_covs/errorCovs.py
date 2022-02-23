# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
import HL_error_covs.arrays as arrays
from HL_error_covs.profiles import ObsProfiles
from HL_error_covs.io_data import IO
from HL_error_covs.utils import Utils
from HL_error_covs.masks import applyMask

# Initialising the classes
ObsProfiles = ObsProfiles()
IO = IO()
Utils = Utils()
applyMask = applyMask()

class HLerrorCovs():
      """ Class of functions to calculate error covariances """
      
      def __init__(self):
          self.fbrmdi = 99999.


      def mp_calc_cov_accum_stats(self, args):
          """ Routine to run on multiprocessors that reads in feedback 
              file and calculates accumulated grid and covariance statistics  
      
              ********* PARAMETERS ********
              1. args: dictionary with members: list of files, bins, 
                       grid_lat, grid_lon, depth_range, obs_type and 
                       source_types

              ********* RETURNS ***********
              1. sum_stats: CovSumStats object containing arrays of
                            summed covariance statistical quantities
                            for each grid box
              2. grid_stats: GridSumStats object containing arrays of
                             within grid summed statistical quantities
                             for each grid box
          """
          # initalises data array classes
          nlat = len(args["grid_lat"])
          nlon = len(args["grid_lon"])
          nbin = len(args["bins"])
          depth_range = args["depth_range"]
          source_types = args["source_types"]
          cov_stats = arrays.CovSumStats((nlat, nlon, nbin))
          grid_stats = arrays.GridSumStats((nlat, nlon))

          # loop over files and calculate the accumulated stats
          for infile in args["list_of_files"]:
              if depth_range:
                  print("Processing file: {} | Depths: {} to {}".format(infile, 
                                     str(depth_range[0]), str(depth_range[1])))
              else:
                  print("Processing file: {} | Depths: Surface variable".format(infile))
      
              # Read fdbk variables
              fdbk_var_array, depths = IO.ncread_fdbk_vars(infile, 
                          args['obs_type'], args['source_types'], args['qc_val'])

              # For profiles pick a single observation at each depth range/latitude/longitude
              # This is done to avoid profiles being correlated with themselves
              if depth_range:
                 fdbk_var_array = ObsProfiles.random_subsample_profiles(depths, 
                                                   depth_range, fdbk_var_array)
                                          
              if len(fdbk_var_array.lats)>0:
                 if (np.min(args["grid_lon"]) >= 0.):
                    fdbk_var_array.lons[np.where(fdbk_var_array.lons < 0.)] = fdbk_var_array.lons + 360.
                           
                 # Need to squash obs and model arrays to be 1D array
                 fdbk_var_array.mod_vals = fdbk_var_array.mod_vals.flatten()
                 fdbk_var_array.obs_vals = fdbk_var_array.obs_vals.flatten()
              
                 # Update stats with summed quantities for each call
                 cov_stats, grid_stats = self.calc_cov_stats(args["grid_lat"], args["grid_lon"], 
                                            args["bins"], fdbk_var_array, cov_stats, grid_stats)
                      
          return cov_stats, grid_stats


      def calc_cov_stats(self, grid_lat, grid_lon, bins, fdbk_var_array, cov_stats, grid_stats):
          """ Routine that acumulates statistics required for covariance 
              calculation from feedback file variables

          ****** PARAMETERS *******
          1. grid_lat, grid_lon: array of latitude and longitude
          2. bins: tuple defining the upper boundary of the bins of separation distance
          3. fdbk_var_array: Object containing the arrays describing the innovations
          4. cov_stats: object containing arrays of summed covariance statistical
                        quantities for each grid box
          5. grid_stats: object containing arrays of summed grid statistical 
                         quantities for each grid box
          
          ******* RETURNS *********
          1. cov_stats: updated object containing arrays of summed covariance
                        statistical quantities for each grid box
          2. grid_stats: updated object containing arrays of summed grid 
                         statistical quantities for each grid box
          """  
          nlat=len(grid_lat)
          nlon=len(grid_lon)
          nbin=len(bins)
          
          lmask = applyMask.mask_bad_fdbk_data(fdbk_var_array)
          num_obs = np.count_nonzero(lmask)
          dlat = np.abs(grid_lat[0] - grid_lat[1])
          dlon = np.abs(grid_lon[0] - grid_lon[1])

          if(num_obs>0):
              for ilat in range(0, nlat):
                  for ilon in range(0, nlon):
                      # First find out which innovations are within seperation distance of the largest binsize  
                      # Add on the distance from the centre to the grid corners
                      grid_dist = Utils.dstance(grid_lon[ilon], grid_lat[ilat], 
                                                np.array([grid_lon[ilon] + dlon/2.]), 
                                                np.array([grid_lat[ilat] + dlat/2.]))
                      obs_dist = Utils.dstance(grid_lon[ilon], grid_lat[ilat], 
                                               fdbk_var_array.lons, fdbk_var_array.lats)
                      lmask1 = np.logical_and(lmask, obs_dist <= bins[nbin-1] + grid_dist[0])
      
                      # Which observations are within the grid square being considered
                      lmask2 = applyMask.mask_obs_in_grid(lmask, ilon, ilat, grid_lon, 
                                                          grid_lat, fdbk_var_array)
                      num_obs_in_grid = np.count_nonzero(lmask2)

                      # calculate grid stats and accumulate them
                      grid_stats = self.calc_grid_stats(ilat, ilon, num_obs_in_grid, 
                                                 fdbk_var_array, lmask2, grid_stats)  

                      # Need at least 2 innovations in a grid box to carry out the 
                      # cov stats calculation 
                      if (num_obs_in_grid >= 2):
                          
                          # Remove those innovations that are in the grid box from ErrsReq to avoid double counting
                          lmask3 = np.logical_and(np.logical_not(lmask2), lmask1)
                          num_independent_obs_req = np.count_nonzero(lmask3)
                           
                          # Set up arrays of innovation using the masks  
                          errs_req = fdbk_var_array.obs_vals[lmask3] - fdbk_var_array.mod_vals[lmask3]
                          lats_req = fdbk_var_array.lats[lmask3]
                          lons_req = fdbk_var_array.lons[lmask3]
      
                          errs_in_grid = fdbk_var_array.obs_vals[lmask2] - fdbk_var_array.mod_vals[lmask2]
                          lats_in_grid = fdbk_var_array.lats[lmask2]
                          lons_in_grid = fdbk_var_array.lons[lmask2] 
      
                          # Pair up innovations in the grid with those outside of it
                          cov_pairs = self.extra_grid_cov_pairs(num_independent_obs_req, num_obs_in_grid, 
                                                           errs_in_grid, errs_req, lons_in_grid, lats_in_grid, 
                                                           lons_req, lats_req)
                          
                          # Next pair-up innovations within the grid box, making sure not to double count 
                          cov_pairs2 = self.intra_grid_cov_pairs(num_obs_in_grid, errs_in_grid, 
                                                                 lons_in_grid, lats_in_grid)
    
                          # Stick these arrays together
                          cov_pairs.append_cov_pair(cov_pairs2)
                          
                          # bin-up the pairs and calculate accumulated cov stats
                          cov_stats = self.accumulate_binned_stats(ilat, ilon, bins, cov_pairs, cov_stats)
         
          return cov_stats, grid_stats


      def extra_grid_cov_pairs(self, num_independent_obs_req, num_obs_in_grid, errs_in_grid, 
                               errs_req, lons_in_grid, lats_in_grid, lons_req, lats_req):
          """ Pairs up innovations in a grid box with given innovations outside of it
          
          ********* PARAMETERS *********
          1. num_independent_obs_req: number of innovations outside of the grid box to be paired
          2. num_obs_in_grid: number of innovations within the grid to be paired
          3. errs_in_grid: array of the innovation values in the grid
          4. errs_req: array of innovations values outside of the grid
          5. lons_in_grid: array of longitudes of the innovations within the grid
          6. lats_in_grid: array of latitudes of the innovations within the grid
          7. lons_req: array of longitudes of the innovations outside of the grid
          8. lats_req: array of latitudes of the innovations outside of the grid 
      
          ********** RETURNS ***********
          1. cov_pairs: object containing arrays of paired statistical quantities
          """
          cov_pairs = arrays.CovPairArrays((num_independent_obs_req*num_obs_in_grid))
          for iob in range(0, num_obs_in_grid):
              start_index = iob*num_independent_obs_req
              end_index = (iob+1)*num_independent_obs_req                       
              cov_pairs.x[start_index:end_index] = errs_in_grid[iob]
              cov_pairs.y[start_index:end_index] = errs_req[:]
              cov_pairs.x_sq[start_index:end_index] = errs_in_grid[iob]*errs_in_grid[iob]
              cov_pairs.y_sq[start_index:end_index] = errs_req[:]*errs_req[:]
              cov_pairs.xy[start_index:end_index] = errs_in_grid[iob]*errs_req[:]        
              cov_pairs.sep_dist[start_index:end_index] = Utils.dstance(lons_in_grid[iob], 
                                              lats_in_grid[iob], lons_req[:], lats_req[:])
          return cov_pairs 


      def intra_grid_cov_pairs(self, num_obs_in_grid, errs_in_grid, lons_in_grid, lats_in_grid):
          """ Pairs up innovations in a grid box with other innovations in the grid
              Double counting is avoided.
          
          ********* PARAMETERS *********
          1. num_obs_in_grid: number of innovations within the grid to be paired
          2. errs_in_grid: array of the innovation values in the grid
          3. lons_in_grid: array of longitudes of the innovations within the grid
          4. lats_in_grid: array of latitudes of the innovations within the grid
      
          ********** RETURNS ***********
          1. cov_pairs: object containing arrays of paired statistical quantities
          """
          grid_cov_elements = (((num_obs_in_grid-1)*(num_obs_in_grid))//2)
          cov_pairs = arrays.CovPairArrays((grid_cov_elements))
          for iob in range(0, num_obs_in_grid-1):
              start_index = ((num_obs_in_grid-1)*(iob)-((iob)*(iob-1))//2)
              end_index = ((num_obs_in_grid-1)*(iob+1)-((iob+1)*(iob))//2)
              cov_pairs.x[start_index:end_index] = errs_in_grid[iob]
              cov_pairs.y[start_index:end_index] = errs_in_grid[iob+1:]
              cov_pairs.x_sq[start_index:end_index] = errs_in_grid[iob]*errs_in_grid[iob]
              cov_pairs.y_sq[start_index:end_index] = errs_in_grid[iob+1:]*errs_in_grid[iob+1:]
              cov_pairs.xy[start_index:end_index] = errs_in_grid[iob]*errs_in_grid[iob+1:]        
              cov_pairs.sep_dist[start_index:end_index] = Utils.dstance(lons_in_grid[iob], lats_in_grid[iob], 
                                                                  lons_in_grid[iob+1:], lats_in_grid[iob+1:])
          return cov_pairs

                              
      def accumulate_binned_stats(self, ilat, ilon, bins, cov_pairs, sum_stats): 
          """ Accumulate the statistics over discrete separation distances between
              the pairs of innovations 
      
          *********** PARAMETERS ***********
          1. ilat, ilon: index of the grid box location    
          2. bins: tuple defining the upper boundary of the bins of separation distance
          3. cov_pairs: object containing arrays of paired statistical quantities
          4. sum_stats: object containing arrays of summed covariance
                        statistical quantities

          *********** RETURNS ************* 
          1. sum_stats: updated object containing arrays of summed covariance
                        statistical quantities
          """
          # Accumulate stats over bins
          for n in range(0, len(bins)):
              lmask = False
              if(n==0):
                lmask = np.logical_and(cov_pairs.sep_dist >= 0., cov_pairs.sep_dist <= bins[n])
              else:
                lmask = np.logical_and(cov_pairs.sep_dist > bins[n-1], cov_pairs.sep_dist <= bins[n])
              num_pairs = np.count_nonzero(lmask)
              if(num_pairs>=2):
                  sum_stats.num_pairs_in_cov[ilat,ilon,n] = sum_stats.num_pairs_in_cov[ilat,ilon,n] + num_pairs
                  sum_stats.sum_x[ilat,ilon,n] = sum_stats.sum_x[ilat,ilon,n] + np.sum(cov_pairs.x[lmask])
                  sum_stats.sum_y[ilat,ilon,n] = sum_stats.sum_y[ilat,ilon,n] + np.sum(cov_pairs.y[lmask])
                  sum_stats.sum_x_sq[ilat,ilon,n] = sum_stats.sum_x_sq[ilat,ilon,n] + np.sum(cov_pairs.x_sq[lmask])
                  sum_stats.sum_y_sq[ilat,ilon,n] = sum_stats.sum_y_sq[ilat,ilon,n] + np.sum(cov_pairs.y_sq[lmask])
                  sum_stats.sum_xy[ilat,ilon,n] = sum_stats.sum_xy[ilat,ilon,n] + np.sum(cov_pairs.xy[lmask])
          return sum_stats


      def calc_grid_stats(self, ilat, ilon, num_obs, fdbk_var_array, lmask, grid_stats):
          """ Function that accumulates within grid statistics from feedback file variables  
      
          ********* PARAMETERS *********
          1. ilat, ilon: index of the grid box location
          2. num_obs: number of observations within the grid box
          3. fdbk_var_array: feedback variables
          4. lmask: logical mask array    
          5. grid_stats: object containing arrays of within grid summed 
                         statistical quantities for each grid box
     
          ********* RETURNS ************
          1. grid_stats: updated object containing arrays of within grid summed 
                         statistical quantities for each grid box
          """
          if(num_obs>0): 
            errs_req = fdbk_var_array.obs_vals[lmask] - fdbk_var_array.mod_vals[lmask]
            grid_stats.num_obs_in_grid[ilat,ilon] = grid_stats.num_obs_in_grid[ilat,ilon] + num_obs
            grid_stats.grid_sum_obs_std[ilat,ilon] = grid_stats.grid_sum_obs_std[ilat,ilon] + \
                                                     np.sum(fdbk_var_array.obs_std[lmask])
            grid_stats.grid_sum[ilat,ilon] = grid_stats.grid_sum[ilat,ilon] + np.sum(errs_req)
            grid_stats.grid_sum_sq[ilat,ilon] = grid_stats.grid_sum_sq[ilat,ilon] + np.sum(errs_req*errs_req)
          return grid_stats
      
     
      def calc_err_covs(self, sum_stats, grid_stats, nbin, nlat, nlon):
          """ Calculate error covariances 
               
              ********* PARAMETERS ********
              1. sum_stats: accumulated covariance statistics
              2. grid_stats: accumulated grid statistics
              2. nbin: number of bins
              3. nlat: number of latitudes
              4. nlon: number of longitudes

              ********* RETURNS ********
              1. cov_xy: covariance for each bin of separation distance and at each grid
              2. corr_xy: correlation for each bin of separation distance and at each grid 
              3. grid_mean: grid mean binned error 
              4. grid_var: grid variance 
              5. num_obs_in_grid: number of obs within squared grid
              6. num_pairs_in_cov: number of covariance pairs for each bin at each grid
              7. stderr_xy: Standard error of the covariance estimate
          """
              
          # Variance of x
          var_x = sum_stats.sum_x_sq/sum_stats.num_pairs_in_cov - \
                 (sum_stats.sum_x/sum_stats.num_pairs_in_cov)* \
                 (sum_stats.sum_x/sum_stats.num_pairs_in_cov)
 
          # Variance of y
          var_y = sum_stats.sum_y_sq/sum_stats.num_pairs_in_cov - \
                 (sum_stats.sum_y/sum_stats.num_pairs_in_cov)* \
                 (sum_stats.sum_y/sum_stats.num_pairs_in_cov)
 
          cov_xy = np.ma.zeros((nlat, nlon, nbin))
          corr_xy = np.ma.zeros((nlat, nlon, nbin))

          # The std error on the variance is var * sqrt(2/(N-1)) for a Guassian
          # Using that (see https://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf)
          stderr_xy = np.ma.zeros((nlat, nlon, nbin))
 
          # Covariance and correlation for each bin
          for n in range(0, nbin):
              cov_xy[:,:,n] = sum_stats.sum_xy[:,:,n]/sum_stats.num_pairs_in_cov[:,:,n] - \
                             (sum_stats.sum_x[:,:,n]/sum_stats.num_pairs_in_cov[:,:,n])* \
                             (sum_stats.sum_y[:,:,n]/sum_stats.num_pairs_in_cov[:,:,n])
              corr_xy[:,:,n] = cov_xy[:,:,n]/np.ma.sqrt(var_x[:,:,n]*var_y[:,:,n])
              stderr_xy[:,:,n] =
 
          # Squared grid mean binned error and variance
          grid_mean = grid_stats.grid_sum/grid_stats.num_obs_in_grid
          grid_var = grid_stats.grid_sum_sq/grid_stats.num_obs_in_grid - \
                    (grid_mean)*(grid_mean)
          numpairscov = sum_stats.num_pairs_in_cov
          numobsingrid = grid_stats.num_obs_in_grid
 
          # Grid mean obs measurement error
          grid_mean_obstd = grid_stats.grid_sum_obs_std/grid_stats.num_obs_in_grid

          return cov_xy, corr_xy, grid_mean, grid_var, numobsingrid, numpairscov, grid_mean_obstd
