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
import modules.arrays as arrays
from modules.io_data import IO
from modules.utils import Utils
from HL_error_covs.errorCovs import HLerrorCovs
from modules.masks import applyMask

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
    outfile = IO.nc_define_dimensions(outfilename,
                                      ['latitude', 'longitude', 'bins', 'depth'],
                                      [nlat, nlon, nbin, ndep])
    
    # Write dimension variables
    IO.ncwrite_variables(outfile, ['latitude'], ['f'], ('latitude'), [grid_lat])
    IO.ncwrite_variables(outfile, ['longitude'], ['f'], ('longitude'), [grid_lon])
    IO.ncwrite_variables(outfile, ['depth'], ['f'], ('depth'), [depths])
    IO.ncwrite_variables(outfile, ['bins'], ['f'], ('bins'), [bins])

    # Create 4D netcdf variables (accumulated stats)
    IO.ncwrite_variables(outfile, ['SumX', 'SumY', 'SumXSq', 'SumYSq', 'SumXY', 'NumObsCov'],
                         ['f', 'f', 'f', 'f', 'f', 'i'], ('depth', 'latitude', 'longitude', 'bins'))

    # Create 3D netcdf variables (accumulated stats)
    IO.ncwrite_variables(outfile, ['GridSum', 'GridSumSq', 'GridSumObsStd', 'GridNumObs'],
                         ['f', 'f', 'f', 'i'], ('depth', 'latitude', 'longitude'))

    print("MESSAGE: {} nprocs to process {} feedback files".format(nproc,
                                               len(list_of_fdbackfiles)))

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
  
        # Write 4D accumulated stats to output file
        IO.ncwrite_variables(outfile, ['SumX', 'SumY', 'SumXSq', 'SumYSq', 'SumXY', 'NumObsCov'], [], [], \
                             vardata=[sum_stats.sum_x, sum_stats.sum_y, sum_stats.sum_x_sq, \
                              sum_stats.sum_y_sq, sum_stats.sum_xy, sum_stats.num_pairs_in_cov], \
                             create_vars=False, dep_lev=dep_lev)

        # Write 3D accumulated stats to output file
        IO.ncwrite_variables(outfile, ['GridSum', 'GridSumSq', 'GridNumObs', 'GridSumObsStd'], [], [], \
                             vardata=[grid_stats.grid_sum, grid_stats.grid_sum_sq, \
                             grid_stats.num_obs_in_grid, grid_stats.grid_sum_obs_std], \
                             create_vars=False, dep_lev=dep_lev)

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

    ncdata = IO.ncread_variables(list_of_files[0], ['latitude', 'longitude', 'bins', 'depth'])
    grid_lat = ncdata[0]
    grid_lon = ncdata[1]
    bins = ncdata[2]
    depths = ncdata[3]
    nlat = len(grid_lat)
    nlon = len(grid_lon)
    nbin = len(bins)
    ndep = len(depths)

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename,
                                      ['latitude', 'longitude', 'bins', 'depth'],
                                      [nlat, nlon, nbin, ndep])
    
    # Write dimension variables
    IO.ncwrite_variables(outfile, ['latitude'], ['f'], ('latitude'), [grid_lat])
    IO.ncwrite_variables(outfile, ['longitude'], ['f'], ('longitude'), [grid_lon])
    IO.ncwrite_variables(outfile, ['depth'], ['f'], ('depth'), [depths])
    IO.ncwrite_variables(outfile, ['bins'], ['f'], ('bins'), [bins])

    # Add netcdf 4D variables (final ErrorCovs stats)
    IO.ncwrite_variables(outfile, ['NumObsCov', 'Covariance', 'Correlation'],
                         ['i', 'f', 'f'], ('depth', 'latitude', 'longitude', 'bins'))

    # Add netcdf 3D variables (final ErrorCovs stats)
    IO.ncwrite_variables(outfile, ['GridVariance', 'GridNumObs', \
                         'GridMeanBinnedError', 'GridMeanObsStd'], ['f', 'i', 'f', 'f'],
                         ('depth', 'latitude', 'longitude'))

    for dep_lev in range(0, ndep):
        print("MESSAGE: Calculating error covariance for level: " + str(depths[dep_lev]) + " m")
        final_cov_stats = arrays.CovSumStats((nlat, nlon, nbin))
        final_grid_stats = arrays.GridSumStats((nlat, nlon))
        for f in list_of_files:
            print("MESSAGE: Reading file {}".format(f))

            # Read cov_stats variables from netcdf
            ncdata = IO.ncread_variables(f, ['SumX', 'SumY', 'SumXSq', 'SumYSq', 'SumXY', \
                                         'NumObsCov'], dep_lev=dep_lev)

            cov_stats = arrays.CovSumStats((nlat, nlon, nbin))
            cov_stats.sum_x = ncdata[0]
            cov_stats.sum_y = ncdata[1]
            cov_stats.sum_x_sq = ncdata[2]
            cov_stats.sum_y_sq = ncdata[3]
            cov_stats.sum_xy = ncdata[4]
            cov_stats.num_pairs_in_cov = ncdata[5]

            # Read grid_stats variables from netcdf
            ncdata = IO.ncread_variables(f, ['GridSum', 'GridSumSq', 'GridNumObs', \
                                         'GridSumObsStd'], dep_lev=dep_lev)

            grid_stats = arrays.GridSumStats((nlat, nlon))
            grid_stats.grid_sum = ncdata[0]
            grid_stats.grid_sum_sq = ncdata[1]
            grid_stats.num_obs_in_grid = ncdata[2]
            grid_stats.grid_sum_obs_std = ncdata[3]

            final_cov_stats += cov_stats
            final_grid_stats += grid_stats

        # Calculate correlation and covariance
        cov_xy, corr_xy, grid_mean, grid_var, numobsgrid, \
        numpairscov, grid_mean_obstd = HLerrorCovs.calc_err_covs(final_cov_stats,
                                              final_grid_stats, nbin, nlat, nlon)

        # Mask output data
        applyMask.mask_output_cov_data(grid_mean, grid_var, grid_mean_obstd,
                                       numobsgrid, numpairscov, cov_xy, corr_xy)

        # Write 4D variables to output file
        IO.ncwrite_variables(outfile, ['NumObsCov', 'Covariance', 'Correlation'], [], [], \
                             vardata=[numpairscov, cov_xy, corr_xy], create_vars=False,
                             dep_lev=dep_lev)

        # Write 3D variables to output file
        IO.ncwrite_variables(outfile, ['GridVariance', 'GridNumObs', 'GridMeanBinnedError', \
                             'GridMeanObsStd'], [], [], vardata=[grid_var, numobsgrid, \
                              grid_mean, grid_mean_obstd], create_vars=False, dep_lev=dep_lev)

    outfile.close()
