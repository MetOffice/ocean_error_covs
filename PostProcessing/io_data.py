# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
from netCDF4 import Dataset

class IO():
    """ Class of IO functions """

    def __init__(self):
        self.MASK_VAL = -1e10


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
            nc = Dataset(infile)
        except:
            raise ValueError("[ERROR] CANNOT OPEN NETCDF FILE")

        grid_lat = nc.variables['latitude'][:]
        grid_lon = nc.variables['longitude'][:]
        depths = nc.variables['depth'][:]
        bins = nc.variables['bins'][:]
        nc.close()
        return grid_lat, grid_lon, depths, bins


    def ncread_errorcovs(self, infile, dep_lev):
        """ Read netcdf with error covariance stats and return
            objects containing variables

        ****** PARAMETERS *****
        1. infile: file to read
        2. dep_lev: depth level to extract variables

        ******* RETURNS *******
        1. var: variance of innovations within each grid box
        2. cors: correlation of innovations within bins of separation distance
        3. gridnumobs: number of obs within each grid box
        """
        # Open netcdf file
        try:
            nc = Dataset(infile)
        except:
            raise ValueError("[ERROR] CANNOT OPEN NETCDF FILE")

        numobsvar = nc.variables["GridNumObs"][dep_lev,:,:]
        var = nc.variables["GridVariance"][dep_lev,:,:]
        cors = nc.variables["Correlation"][dep_lev,:,:,:]

        # account for precision errors by forcing minimum and 
        # maximum correlation to [-1.0, 1.0]
        cors[cors>1.] = 1.
        cors[cors<-1.] = -1.
        nc.close()
        return var, cors, numobsvar


    def nc_define_dimensions(self, outfilename, nlat, nlon, ndep):
        """ Define output netcdf file 

        ****** PARAMETERS *****
        1. outfilename: path with filename
        2. nlat: number of latitudes
        3. nlon: number of longitudes
        4. ndep: number of depth levels

        ******* RETURNS *******
        1. outfile: object containing netcdf variables and attributes. 
        """ 
  
        # Open netcdf file object 
        outfile = Dataset(outfilename, mode='w', clobber=True)

        # define dimensions
        outfile.createDimension('latitude', size=nlat)
        outfile.createDimension('longitude', size=nlon)
        outfile.createDimension('depth', size=ndep)
        
        # define variables
        outfile.createVariable("latitude",'f', dimensions=("latitude"))
        outfile.createVariable("longitude",'f', dimensions=("longitude"))
        outfile.createVariable("depth", 'f', dimensions=("depth"))
        return outfile


    def nc_define_vars(self, outfile, varstring, vartype, dimensions):
        """ Define variables to be written in output netcdf file

        ****** PARAMETERS *****
        1. outfile: object containing netcdf variables and attributes
        2. varstring: name of the variable to appear in the netcdf file
        3. vartype: variable type 
        4. dimensions: tuple with name of each dimension
        """ 
        outfile.createVariable(varstring, vartype, dimensions=dimensions)


    def ncwrite_dimension_variables(self, outfile, grid_lat, grid_lon, depths):
        """ Write variables that are part of the netcdf dimensions 
            (e.g. lat, lon and depths)

        ****** PARAMETERS *****
        1. outfile: handle to the output file
        2. grid_lat: latitudes
        3. grid_lon: longitudes
        4. depth_levels: depth levels (for profile only) 
        """ 
        outfile.variables["latitude"][:] = grid_lat
        outfile.variables["longitude"][:] = grid_lon
        outfile.variables["depth"][:] = depths


    def ncwrite_output(self, outfile, func, chi_grid, obs_err, params, lev):
        """ Write variables of netcdf file per depth level
        
        ******* PARAMETERS *******
        1. outfile: Netcdf object 
        2. func: instance of fitting function
        3. chi_grid: gridded chi squared criterion
        4. obs_err: gridded obs error 
        5. params: results of fitting function
        6. lev: depth level 
        """ 

        # masking variables
        msk = np.logical_or(chi_grid == self.MASK_VAL, np.isnan(chi_grid))
        chi_grid.mask = msk
        obs_err.mask  = msk
        for param in range(0, len(params)):
            params[param].mask = msk   
        
        # write data per depth level
        outfile.variables["Chi_sq"][lev,:,:] = chi_grid
        outfile.variables["obs_err"][lev,:,:] = obs_err
        for param in range(0, len(params)):
            outfile.variables[func.param_names()[param]][lev,:,:] = params[param][:,:]

