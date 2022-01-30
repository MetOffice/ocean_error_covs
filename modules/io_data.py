# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
from netCDF4 import Dataset
from modules.masks import applyMask
import modules.arrays as arrays

# Initialising the classes
applyMask = applyMask()

class IO():
    """ Class of IO functions """

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
                          vardata=[], dep_lev=None, create_vars=True):
        """ Create and write netcdf variables

        ****** PARAMETERS *****
        1. outfile: netcdf object
        2. varnames: names of the variables
        3. vartypes: types of the variables (float, integer, etc)
        4. vardims: tuple with the dimension names
        5. vardata: list containing all the arrays to be written
        6. dep_lev: slice of depth level to be written
        7. create_vars: True if step to create vars should be done
        """ 

        # create dimension and write variables
        for idx, var in enumerate(varnames):
            if create_vars:
               outfile.createVariable(var, vartypes[idx], dimensions=vardims)
            if len(vardata) > 0:
               if dep_lev is None:
                  outfile.variables[var][:] = vardata[idx]
               else:
                  outfile.variables[var][dep_lev,:] = vardata[idx]

