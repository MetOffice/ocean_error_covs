# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
################## Python packages ###########################
import os
# Set os.environment to have NUM_THREADS="1" to avoid numpy using 
# by default all processors and competing with parallel processes of 
# error covariance calculation.
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from multiprocessing import Pool
################## Code modules ##############################
import HL_error_covs.arrays as arrays
from HL_error_covs.io_data import IO
from HL_error_covs.utils import Utils
from HL_error_covs.errorCovs import HLerrorCovs
from HL_error_covs.masks import applyMask

# Initialising the classes
IO = IO()
Utils = Utils()
HLerrorCovs = HLerrorCovs()
applyMask = applyMask()

def HL_cov_accum_stats(list_of_fdbackfiles, obs_type="SST",
                       outfilename="accum_stats.nc", grid_def=[[-90,90,2.],[-180,180,2.]],
                       bins=np.arange(50,1000.,50.), depth_boundaries=[],
                       source_types=[], qc_val=[], nproc=1):

    """ Top-level routine that calculates the H+L accumulated statistics from a 
        list of feedback files using multiprocessing and produces an output file
        which is further used to calculate the error covariances.

    *************** PARAMETERS *******************
    1. list_of_fdbackfiles: list of feedback files
    2. outfilename: name of file to contain the output statistics
    3. grid_def: Definition of pre-sorting grid; should be a list-of-lists of the form
                 [[min lat, max lat, delta lat], [min lat, max lat, delta lat]]
    4. obs_type: observation type to process
    5. bins: list defining the upper boundary of the bins of separation distance
             (in km) to be used for correlation calculation
    6. depth_boundaries: depth level boundaries (only used for processing 
             profile observations)
    7. source_types: observation id types to process 
                     (default: [], which means process all)
    8. qc_val: QC values to filter obs
                     (default: [], which means do not filter any obs)
    9. nproc: number of processors (default is 1)
    """ 
    # Check list of feedack files
    list_of_fdbackfiles = Utils.check_files(list_of_fdbackfiles)
    if list_of_fdbackfiles is None:
       raise ValueError("[ERROR] FEEDBACK FILES NOT FOUND")

    # define number of bins
    nbin = len(bins)
    bins = np.array(bins)

    # define depth variables 
    ndep = 1
    depths = 0
    if any(depth_boundaries):
       depth_boundaries = np.array(depth_boundaries)
       ndep = len(depth_boundaries)-1
       depths = 0.5 * (depth_boundaries[:-1] + depth_boundaries[1:])
    
    # create variables for pre sorting grid
    grid_lat = np.arange(grid_def[0][0],grid_def[0][1],grid_def[0][2]) + grid_def[0][2]/2.
    grid_lon = np.arange(grid_def[1][0],grid_def[1][1],grid_def[1][2]) + grid_def[1][2]/2.
    nlat = len(grid_lat)
    nlon = len(grid_lon)

    # divide list of files for each processor
    list_per_proc = Utils.divide_files_per_proc(nproc, list_of_fdbackfiles)

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename, nlat, nlon, nbin, ndep)
    
    # Add netcdf variables (accumulated stats)
    IO.nc_define_accum_stats_variables(outfile)

    # Write dimension variables (depth, lat, lon and bins)
    IO.ncwrite_dimension_variables(outfile, grid_lat, grid_lon, depths, bins)

    print("MESSAGE: {} nprocs to process {} feedback files".format(nproc,
                                               len(list_of_fdbackfiles)))
    arg_list=[]
    for dep_lev in range(0, ndep):
        print("MESSAGE: Calculating accumulated stats for level: ", dep_lev)
        
        # set up workers
        workers=Pool(nproc)
        
        if ndep == 1:
            depth_range = []
        else:
            depth_range = [depth_boundaries[dep_lev],depth_boundaries[dep_lev+1]]
            
        arg_list=[]
        for n in range(0, nproc):
            arg_list += [{"list_of_files":list_per_proc[n],
                          "bins": bins*1000,    # NOTE: conversion to meters
                          "grid_lon": grid_lon,
                          "grid_lat": grid_lat,
                          "depth_range": depth_range,
                          "obs_type": obs_type,
                          "source_types": source_types,
                          "qc_val": qc_val}]

        # send tasks off to workers
        work_output = workers.map(HLerrorCovs.mp_calc_cov_accum_stats, arg_list)
        workers.close()
 
        # Accumulate stats over all the processors
        sum_stats = arrays.CovSumStats((nlat, nlon, nbin))
        grid_stats = arrays.GridSumStats((nlat, nlon))
        for p in work_output:
            sum_stats += p[0]
            grid_stats += p[1]
  
        # Write accumulated stats to output file
        IO.ncwrite_accum_stats(outfile, dep_lev, sum_stats, grid_stats)

    outfile.close()


def HL_error_covs(list_of_files, outfilename="corrs.nc"):

    """ Top-level routine that calculates the H+L error covariances from a 
        list of netcdf files containing the accumulated statistics for each
        grid box.

    *************** PARAMETERS *******************
    1. list_of_files: list of netcdf files containing the accumulated statistics
    2. outfilename: name of file to contain the output statistics
    """ 
    # Check list of netcdf files
    list_of_files = Utils.check_files(list_of_files)
    if list_of_files is None:
       raise ValueError("[ERROR] NETCDF FILES NOT FOUND")

    grid_lat, grid_lon, depths, bins = IO.ncread_dimension_variables(list_of_files[0])
    nbin = len(bins)
    ndep = len(depths)
    nlat = len(grid_lat)
    nlon = len(grid_lon)

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename, nlat, nlon, nbin, ndep)
    
    # Add netcdf variables (final ErrorCovs stats)
    IO.nc_define_cov_variables(outfile)

    # Write dimension variables (depth, lat, lon and bins)
    IO.ncwrite_dimension_variables(outfile, grid_lat, grid_lon, 
                                   depths, bins)

    for dep_lev in range(0, ndep):
        print("MESSAGE: Calculating error covariance for level: " + str(depths[dep_lev]) + " m")
        final_cov_stats = arrays.CovSumStats((nlat, nlon, nbin))
        final_grid_stats = arrays.GridSumStats((nlat, nlon))
        for f in list_of_files:
            print("MESSAGE: Reading file {}".format(f))
            cov_stats, grid_stats = IO.ncread_accum_stats(f, nlat, nlon, nbin, dep_lev)
            final_cov_stats += cov_stats
            final_grid_stats += grid_stats

        # Calculate correlation and covariancee
        cov_xy, corr_xy, grid_mean, grid_var, numobsgrid, \
        numpairscov, grid_mean_obstd = HLerrorCovs.calc_err_covs(final_cov_stats,
                                              final_grid_stats, nbin, nlat, nlon)

        # Mask output data
        applyMask.mask_output_cov_data(grid_mean, grid_var, grid_mean_obstd,
                                       numobsgrid, numpairscov, cov_xy, corr_xy)

        # Write error covariances to output file
        IO.ncwrite_covariance(outfile, dep_lev, grid_mean, grid_var, grid_mean_obstd,
                              numobsgrid, numpairscov, cov_xy, corr_xy)

    outfile.close()
