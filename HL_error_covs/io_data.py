# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
from netCDF4 import Dataset
from masks import applyMask
import arrays

# Initialising the classes
applyMask = applyMask()

class IO():
    """ Class of IO functions """

    def __init__(self):
        self.fill_value = 1e10
        self.fill_value_numobs = -99999


    def ncread_fdbk_vars(self, infile, var_type, source_types):
        """ Read netcdf of feedback variables and return
            object containing variables
        
        ****** PARAMETERS *****
        1. infile: netcdf file to be read
        2. var_type: type of observation to be considered
        3. source_types: id of each sensor to be considered

        ******* RETURNS *******
        1. fdbk_var_array: object containing arrays
                           of feedback file variables
        """
        # Open netcdf file
        try:
            fdbkdata = Dataset(infile)
        except:
            raise ValueError("[ERROR] CANNOT OPEN NETCDF FILE")

        fdbk_var_array = arrays.FdbkVarArrays()

        # station type mask
        msk = applyMask.vartype_mask(fdbkdata, source_types)
        if not(np.any(msk)):
           print("[WARNING] OBSERVATION IDs NOT FOUND IN {}".format(infile))
           return fdbk_var_array, None
        try:
            fdbk_var_array.obs_vals = fdbkdata.variables[var_type+"_OBS"][msk,:]
        except:
            raise ValueError("[ERROR] " + var_type + "NOT FOUND IN FEEDBACK FILE " + infile)
        fdbk_var_array.lats = fdbkdata.variables["LATITUDE"][msk]
        fdbk_var_array.lons = fdbkdata.variables["LONGITUDE"][msk]
        fdbk_var_array.obs_qc = fdbkdata.variables[var_type + "_QC"][msk]
        fdbk_var_array.mod_vals = fdbkdata.variables[var_type + "_Hx"][msk,:]
        depths = fdbkdata.variables["DEPTH"][msk,:]
        fdbkdata.close()
        return fdbk_var_array, depths


    def ncread_dimension_variables(self, infile):
        """ Read dimensions of the netcdf file

        ****** PARAMETERS *****
        1. infile: file to read
        
        ******* RETURNS *******
        1. latitude: latitudes of the grid
        2. longitude: longitudes of the grid
        3. depths: depth levels
        4. bins: bins of separation distance
        """
        # Open netcdf file
        try:
            nc = Dataset(infile, mode='r')
        except:
            raise ValueError("[ERROR] CANNOT OPEN NETCDF FILE")

        grid_lat = nc.variables['latitude'][:]
        grid_lon = nc.variables['longitude'][:]
        depths = nc.variables['depth'][:]
        bins = nc.variables['bins'][:]
        nc.close()
        return grid_lat, grid_lon, depths, bins


    def ncread_accum_stats(self, infile, nlat, nlon, nbin, dep_lev):
        """ Read netcdf of accumulated statistics and return
            objects containing variables

        ****** PARAMETERS *****
        1. infile: file to read
        2. nlat: number of latitudes
        3. nlon: number of longitudes
        4. nbin: number of bins of separation distance
        5. dep_level: depth level to extract accumulated statistics

        ******* RETURNS *******
        1. cov_stats: object containing arrays of summed covariance statistical
                      quantities for each grid box
        2. grid_stats: object containing arrays of summed grid statistical 
                       quantities for each grid box
        """
        # Open netcdf file
        try:
            accumstat = Dataset(infile)
        except:
            raise ValueError("[ERROR] CANNOT OPEN NETCDF FILE")

        # Initialise arrays
        cov_stats = arrays.CovSumStats((nlat, nlon, nbin))
        grid_stats = arrays.GridSumStats((nlat, nlon))

        cov_stats.sum_x = accumstat.variables["SumX"][dep_lev,:,:,:]
        cov_stats.sum_y = accumstat.variables["SumY"][dep_lev,:,:,:]
        cov_stats.sum_x_sq = accumstat.variables["SumXSq"][dep_lev,:,:,:]
        cov_stats.sum_y_sq = accumstat.variables["SumYSq"][dep_lev,:,:,:]
        cov_stats.sum_xy = accumstat.variables["SumXY"][dep_lev,:,:,:]
        cov_stats.num_pairs_in_cov = accumstat.variables["NumObsCov"][dep_lev,:,:,:]

        grid_stats.grid_sum = accumstat.variables["GridSum"][dep_lev,:,:]
        grid_stats.grid_sum_sq = accumstat.variables["GridSumSq"][dep_lev,:,:]
        grid_stats.num_obs_in_grid = accumstat.variables["GridNumObs"][dep_lev,:,:]
        accumstat.close()
        return cov_stats, grid_stats


    def nc_define_dimensions(self, outfilename, nlat, nlon, nbin, ndep):
        """ Define output netcdf file 

        ****** PARAMETERS *****
        1. outfilename: path with filename
        2. nlat: number of latitudes
        3. nlon: number of longitudes
        4. nbin: number of correlation bins
        5. ndep: number of depth levels

        ******* RETURNS *******
        1. outfile: object containing netcdf variables and attributes. 
        """ 
  
        # Open netcdf file object 
        outfile = Dataset(outfilename, mode='w', clobber=True)

        # define dimensions
        outfile.createDimension('latitude', size=nlat)
        outfile.createDimension('longitude', size=nlon)
        outfile.createDimension('bins', size=nbin)
        outfile.createDimension('depth', size=ndep)
        
        # define variables
        outfile.createVariable("latitude",'f', dimensions=("latitude"))
        outfile.createVariable("longitude",'f', dimensions=("longitude"))
        outfile.createVariable("depth", 'f', dimensions=("depth"))
        outfile.createVariable("bins", 'f', dimensions=("bins"))
        return outfile


    def nc_define_accum_stats_variables(self, outfile):
        """ Define accumulated stats to be written in output netcdf file

        ****** PARAMETERS *****
        1. outfile: object containing netcdf variables and attributes
        """ 
        outfile.createVariable("SumX", 'f', dimensions=("depth", "latitude", "longitude", "bins"))
        outfile.createVariable("SumY", 'f', dimensions=("depth", "latitude", "longitude", "bins"))
        outfile.createVariable("SumXSq", 'f', dimensions=("depth", "latitude", "longitude", "bins"))
        outfile.createVariable("SumYSq", 'f', dimensions=("depth", "latitude", "longitude", "bins"))
        outfile.createVariable("SumXY", 'f', dimensions=("depth", "latitude", "longitude", "bins"))
        outfile.createVariable("NumObsCov", 'i', dimensions=("depth", "latitude", "longitude", "bins"))
        outfile.createVariable("GridSum", 'f', dimensions=("depth", "latitude", "longitude"))
        outfile.createVariable("GridSumSq", 'f', dimensions=("depth", "latitude", "longitude"))
        outfile.createVariable("GridNumObs", 'i', dimensions=("depth", "latitude", "longitude"))


    def nc_define_cov_variables(self, outfile):
        """ Define covariance variables to be written in output netcdf file

        ****** PARAMETERS *****
        1. outfile: object containing netcdf variables and attributes
        """ 
        outfile.createVariable("GridVariance", 'f', dimensions=("depth", "latitude", "longitude"), 
                               fill_value=self.fill_value)
        outfile.createVariable("NumObsCov", 'i', dimensions=("depth", "latitude", "longitude", "bins"),
                               fill_value=self.fill_value_numobs)
        outfile.createVariable("GridNumObs", 'i', dimensions=("depth", "latitude", "longitude"),
                               fill_value=self.fill_value_numobs)
        outfile.createVariable("GridMeanBinnedError", 'f', dimensions=("depth", "latitude", "longitude"),
                               fill_value=self.fill_value)
        outfile.createVariable("Covariance", 'f', dimensions=("depth", "latitude", "longitude", "bins"),
                               fill_value=self.fill_value)
        outfile.createVariable("Correlation", 'f', dimensions=("depth","latitude", "longitude", "bins"), 
                                 fill_value=1.e10)


    def ncwrite_dimension_variables(self, outfile, grid_lat, grid_lon, 
                                    depths, bins):
        """ Write variables that are part of the netcdf dimensions 
            (e.g. lat, lon, depths and bins)

        ****** PARAMETERS *****
        1. outfile: handle to the output file
        2. grid_lat: latitudes
        3. grid_lon: longitudes
        4. depth_levels: depth levels (for profile only) 
        5. bins: bins of separation distance
        """ 
        outfile.variables["latitude"][:] = grid_lat
        outfile.variables["longitude"][:] = grid_lon
        outfile.variables["bins"][:] = bins
        outfile.variables["depth"][:] = depths


    def ncwrite_accum_stats(self, outfile, dep_lev, sum_stats, grid_stats):
        """ Write accumulated statistics to netcdf file

        ****** PARAMETERS ******
        1. outfile: handle to the output file
        2. dep_lev: the depth level of the data
        3. sum_stats: object containing summed covariance statistical quantities
        4. grid_stats: object containing summed grid statistical quantities
        """
        outfile.variables["SumX"][dep_lev,:,:,:] = sum_stats.sum_x
        outfile.variables["SumY"][dep_lev,:,:,:] = sum_stats.sum_y
        outfile.variables["SumXSq"][dep_lev,:,:,:] = sum_stats.sum_x_sq
        outfile.variables["SumYSq"][dep_lev,:,:,:] = sum_stats.sum_y_sq
        outfile.variables["SumXY"][dep_lev,:,:,:] = sum_stats.sum_xy
        outfile.variables["NumObsCov"][dep_lev,:,:,:] = sum_stats.num_pairs_in_cov
        outfile.variables["GridSum"][dep_lev,:,:] = grid_stats.grid_sum
        outfile.variables["GridSumSq"][dep_lev,:,:] = grid_stats.grid_sum_sq
        outfile.variables["GridNumObs"][dep_lev,:,:] = grid_stats.num_obs_in_grid


    def ncwrite_covariance(self, outfile, dep_lev, grid_mean, grid_var, 
                          num_obs_in_grid, num_pairs_in_cov, cov_xy, corr_xy):
        """ Write covariance data to output netcdf file 

        ****** PARAMETERS ******
        1. outfile: handle to the output file
        2. dep_lev: the depth level of the data
        3. grid_mean: grid mean binned error 
        4. grid_var: grid variance 
        5. num_obs_in_grid: Number of obs within squared grid
        6. num_pairs_in_cov: Number of covariance pairs for each bin at each grid
        7. cov_xy: Covariance for each bin of separation distance and at each grid
        8. corr_xy: Correlation for each bin of separation distance and at each grid 
        """
        outfile.variables["GridMeanBinnedError"][dep_lev,:,:] = grid_mean
        outfile.variables["GridVariance"][dep_lev,:,:] = grid_var
        outfile.variables["GridNumObs"][dep_lev,:,:] = num_obs_in_grid
        outfile.variables["NumObsCov"][dep_lev,:,:,:] = num_pairs_in_cov
        outfile.variables["Covariance"][dep_lev,:,:,:] = cov_xy
        outfile.variables["Correlation"][dep_lev,:,:,:] = corr_xy
