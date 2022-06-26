# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
import xarray as xr
from datetime import datetime as dt
from netCDF4 import Dataset
from modules.masks import applyMask
import modules.arrays as arrays

# Initialising the classes
applyMask = applyMask()

class IO_netCDF4():
    """ Class of IO functions based on netCDF4 capabilities """

    def __init__(self):
        self.fill_value = 1e10
        self.fill_value_numobs = -99999


    def ncread_fdbk_vars(self, infile, var_type, source_types, qc_val):
        """ Read netcdf of feedback variables and return
            object containing variables
        
        ****** PARAMETERS *****
        1. infile: netcdf file to be read
        2. var_type: type of observation to be considered
        3. source_types: id of each sensor to be considered
        4. qc_val: qc values to filter observations

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
            fdbk_var_array.obs_qc = fdbkdata.variables[var_type + "_QC"][:]
        except:
            raise ValueError("[ERROR] " + var_type + "NOT FOUND IN FEEDBACK FILE " + infile)

        # Filtering observations based on QC flags
        if qc_val:
           for qc in qc_val:
               msk[fdbk_var_array.obs_qc == qc] = False
        fdbk_var_array.obs_qc = fdbk_var_array.obs_qc[msk]

        fdbk_var_array.obs_vals = fdbkdata.variables[var_type + "_OBS"][msk,:]
        try:
            fdbk_var_array.obs_std = (
                fdbkdata.variables[var_type + "_STD"][msk,:])
        except KeyError:
            print("No STD in file, setting fdbk_var_array.obs_std to zero")
            fdbk_var_array.obs_std = np.zeros(fdbk_var_array.obs_vals.shape)

        fdbk_var_array.lats = fdbkdata.variables["LATITUDE"][msk]
        fdbk_var_array.lons = fdbkdata.variables["LONGITUDE"][msk]
        # Handle longitude ambiguity
        fdbk_var_array.lons[fdbk_var_array.lons>180] -= 360
        fdbk_var_array.mod_vals = fdbkdata.variables[var_type + "_Hx"][msk,:]
        depths = fdbkdata.variables["DEPTH"][msk,:]
        fdbkdata.close()
        return fdbk_var_array, depths


    def ncread_variables(self, infile, varnc, dep_lev=None):
        """ Read variables of a netcdf file

        ****** PARAMETERS *****
        1. infile: file to read
        2. varnc: list with varnames to be read
        3. dep_lev: slice of depth level to be read
        
        ******* RETURNS *******
        1. data: data group with all requested variables
        """
        # Open netcdf file
        try:
            nc = Dataset(infile, mode='r')
        except:
            raise ValueError("[ERROR] CANNOT OPEN NETCDF FILE")

        # Read data group based on list of variable names
        data = []
        for var in varnc:
            if dep_lev is None:
                data.append(nc.variables[var][:])
            else:
                data.append(nc.variables[var][dep_lev,:])
        nc.close()

        return data


    def nc_define_dimensions(self, outfilename, vardim_names, vardims):
        """ Define dimensions of output netcdf file

        ****** PARAMETERS *****
        1. outfilename: path with filename
        2. vardim_nanmes: names of the dimension variables
        3. vardims: dimensions of each variable

        ******* RETURNS *******
        1. outfile: object containing netcdf variables and attributes. 
        """ 
  
        # Open netcdf file object 
        outfile = Dataset(outfilename, mode='w', clobber=True)

        # define dimensions
        for idx, var in enumerate(vardim_names):
            outfile.createDimension(var, size=vardims[idx])
        
        return outfile


    def ncwrite_variables(self, outfile, varnames, vartypes, vardims,
                          vardata=[], dep_lev=None, create_vars=True,
                          fill_value=None):
        """ Create and write netcdf variables

        ****** PARAMETERS *****
        1. outfile: netcdf object
        2. varnames: names of the variables
        3. vartypes: types of the variables (float, integer, etc)
        4. vardims: tuple with the dimension names
        5. vardata: list containing all the arrays to be written
        6. dep_lev: slice of depth level to be written
        7. create_vars: True if step to create vars should be done
        8. fill_value: fill value to be considered (otherwise None)
        """ 

        # create dimension and write variables
        for idx, var in enumerate(varnames):
            if create_vars:
               if fill_value is None:
                  outfile.createVariable(var, vartypes[idx], dimensions=vardims)
               else:
                  outfile.createVariable(var, vartypes[idx], dimensions=vardims,
                                         fill_value=fill_value[idx])
            if len(vardata) > 0:
               if dep_lev is None:
                  outfile.variables[var][:] = vardata[idx]
               else:
                  outfile.variables[var][dep_lev,:] = vardata[idx]


class IO_xarray():
    """ Class of IO functions based on xarray capabilities """

    def __init__(self, PARALLEL_IO=True, CHUNKS_INPUT={'time_counter': 1},
                 field_name='votemper', level=None, COMPRESS_OUTPUT=True,
                 depth_names=["deptht", "depthu", "depthv"],
                 time_dims="time_counter", time_var="time_instant",
                 coord_varnames=["nav_lat", "nav_lon", "deptht", "depthu",
                                 "depthv", "time_instant", "deptht_bounds",
                                 "depthu_bounds", "depthv_bounds", "bounds_lon",
                                 "bounds_lat", "area", "time_instant_bounds"]):
        self.PARALLEL_IO = PARALLEL_IO
        self.CHUNKS_INPUT = CHUNKS_INPUT
        self.COMPRESS_OUTPUT = COMPRESS_OUTPUT
        self.field_name = field_name
        self.depth_names = depth_names
        self.coord_varnames = coord_varnames
        self.level = level
        self.time_dims = time_dims
        self.time_var = time_var


    def findvardepth(self, varlist):
        """
        Check from a list of variables which one matches the depth variables

        ****** PARAMETERS *****
        1. varlist: list of netcdf variables

        ******* RETURNS *******
        1. var: target depth variable name
        """
        for var in self.depth_names:
            if var in varlist:
                return var
        return None


    def get_shape_cov(self, data):
        """
        Get the shape of the covariance array

        ****** PARAMETERS *****
        1. data: xarray data
        ******* RETURNS *******
        1. dshape: shape of error covariance array
                   without time dimension
        """
        if self.level is not None:
            data = data.squeeze()
        dshape = data[self.field_name].shape[1::]
        if len(dshape) > 3:
           raise ValueError("Model data has more than 3 dims!")
        return dshape


    def get_input_model_data(self, data, rec_number):
        """
        Get the input data and its validity time

        ****** PARAMETERS *****
        1. data: xarray data
        2. rec_number: time index

        ******* RETURNS *******
        1. rec: model slice based on time record
        2. tim: time record
        """
        if data[self.field_name].dims[0] != self.time_dims:
            raise ValueError("1st dimension of {} is not {} !".format(
                self.field_name, self.time_var))
        rec = data[self.field_name][[rec_number]]
        tim = data[self.time_var][rec_number]
        tim = dt.strptime(np.datetime_as_string(tim, unit='m'),
                          "%Y-%m-%dT%H:%M")
        return rec, tim


    def read_model_files(self, files, concatenate=False, concat_dim=None):
        """
        Read netcdf model file for a particular level from the
        input file(s) using xarray

        ****** PARAMETERS *****
        1. files:   list of files
        2. concatenate: logical for concatenation of files
        3. concat_dim: dimension name to concatenate

        ******* RETURNS *******
        1. An xarray dataset over multiple files for the
           selected variable and level.
        """
        if concatenate:
            input_ds = xr.open_mfdataset(files, concat_dim=concat_dim, combine="nested",
                                         parallel=self.PARALLEL_IO, chunks=self.CHUNKS_INPUT)
        else:
            input_ds = xr.open_mfdataset(files, parallel=self.PARALLEL_IO, chunks=self.CHUNKS_INPUT)

        # Depth variable
        varlist = list(input_ds.variables.keys())
        depthvar = self.findvardepth(varlist)

        if depthvar is None:
            raise ValueError('Depth variable not found in the nc file!')
        if self.level is None:
            data = input_ds.isel()
        else:
            indexers = {depthvar: [self.level]}
            data = input_ds.isel(**indexers)
        return data, depthvar


    def write_forecast_errors(self, data, diff, rec, output_file):
        """
        Write out the model differences to a netcdf file

        ****** PARAMETERS *****
        1. files: xarray data
        2. diff: model differences
        3. rec: time index
        4. output_file: name of netcdf output

        """
        indexers = {self.time_dims: [rec]}
        data = data.isel(**indexers)
        data[self.field_name] = diff

        ncvars = self.coord_varnames
        ncvars.append(self.field_name)
        for var in data.variables.keys():
            if var not in ncvars:
                data = data.drop(var)

        if self.COMPRESS_OUTPUT:
            encoding = {}
            for var in list(data.variables.keys()):
                if var not in self.coord_varnames:
                    encoding.update({var: {"zlib": True,
                                        "complevel": 1}})
        else:
            encoding = None

        data.to_netcdf(output_file, unlimited_dims=[self.time_dims],
                       encoding=encoding)
