# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
################### Python packages ###########################
import os
# Set os.environment to have NUM_THREADS="1" to avoid numpy using 
# by default all processors and competing with parallel processes.
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from multiprocessing import Pool
################## Code modules ##############################
from modules.io_data import IO
from modules.plot import Plots
from modules.masks import applyMask
from modules.posproc import Posproc

# Initialising the classes
IO = IO()
Posproc = Posproc()
Plots = Plots()
applyMask = applyMask()

def HL_fitting_function(infile, outfilename, func_name="MultiGauss", 
                       num_funcs=2, lenscale=(400,40), plot=None, outfig='./figures', 
                       nproc=4, min_num_obs=2, max_iter=100, scalefac=1.0):

    """ Top-level routine that fits a specific function to HL stats covariance file

    ***************** PARAMETERS *******************
    1. infile: name of file containing the HL error covariances
    2. outfilename: name of file to write the results to
    3. func_name: name of function to fit to (options: MultiGauss and MultiGauss_Fixed)
    4. num_funcs: number of functions to use (default: 2)
    5. lenscale: Tuple of pre-defined lengthscales used in MultiGauss_Fixed
                 or in MultiGauss as initial guesses for the lengthscales.
                 Number of tuple members must be equal to number of functions.
    6. plot: positions of (x,y) pairs to plot or None (default: None)
    7. outfig: path to save the figs
    8. nproc: number of processors to use (default: 4)
    9. min_num_obs: minimum number of observations to do calculations
    10. max_iter: max number of iterations
    11. scalefac: factor to scale the variances when they are very small
    """

    # Checking consistency of input parameters
    if(func_name != "MultiGauss" and func_name != "MultiGauss_Fixed"):
       raise ValueError("[ERROR] FUNCTION NOT AVAILABLE")

    if(len(lenscale) != num_funcs):
       raise ValueError("[ERROR] NUMBER OF LENGTHSCALES NOT COMPATIBLE " +
                         "WITH NUMBER OF FUNCTIONS")

    if(not os.path.exists(infile)):
       raise ValueError("[ERROR] INPUT FILE NOT FOUND")
  
    # Read dimension variables from netcdf file
    ncdata = IO.ncread_variables(infile, ['latitude', 'longitude', 'bins', 'depth'])
    lats = ncdata[0]
    lons = ncdata[1]
    bins = ncdata[2]
    depth = ncdata[3]

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename,
                                      ['latitude', 'longitude' , 'depth'],
                                      [len(lats), len(lons), len(depth)])

    # Write dimension variables
    IO.ncwrite_variables(outfile, ['latitude'], ['f'], ('latitude'), vardata=[lats])
    IO.ncwrite_variables(outfile, ['longitude'], ['f'], ('longitude'), vardata=[lons])
    IO.ncwrite_variables(outfile, ['depth'], ['f'], ('depth'), vardata=[depth])

    # Add attributes
    outfile.Function = "Function fitting done using the "+func_name+" function"
        
    # Add variables
    IO.ncwrite_variables(outfile, ['Chi_sq', 'obs_err'],
                         ['f', 'f'], ('depth', 'latitude', 'longitude'))

    # Calculate x positions based on the separation distances
    x_val = Posproc.calc_x_positions(bins)

    for lev in range(0, len(depth)):
        print(f"MESSAGE: Fitting function {func_name} to ErrorCov data: {depth[lev]} m")
        
        # set up workers
        workers = Pool(nproc)
        
        # Reading error covariance variables
        ncdata = IO.ncread_variables(infile, ['GridNumObs', 'GridVariance', \
                                     'Correlation', 'Covariance'], dep_lev=lev)
        numobsvar = ncdata[0]
        var = ncdata[1]
        cors = ncdata[2]
        covs = ncdata[3]

        # account for precision errors by forcing minimum and
        # maximum correlation to [-1.0, 1.0]
        cors[cors>1.] = 1.
        cors[cors<-1.] = -1.

        # Creating list with arguments to run in parallel
        arg_lists = Posproc.create_arg_list(x_val, cors, var, numobsvar, min_num_obs, 
                                            func_name, num_funcs, lenscale, max_iter,
                                            scalefac)
        
        # Get workers to do parallel calculations
        results = workers.map(Posproc.fitter, arg_lists)
        workers.close()

        # Unravel results into output grids
        params, obs_err, chi_grid = Posproc.results_to_grid(results, len(lats), len(lons))

        # Remove scale factor if needed
        obs_err = obs_err/scalefac
        for param in range(0, len(params)):
            if "Magnitude" in arg_lists[0]["func"].param_names()[param]:
               params[param] = params[param]/scalefac

        # Plot some results if requested
        if plot != None:
           print(f"MESSAGE: Plotting results - data versus fitting: {depth[lev]} m")
           Plots.plot_data_vs_fitting(outfig, plot, x_val, cors, var, obs_err, lats, lons,
                                       depth[lev], params, func_name, num_funcs, lenscale)

        print(f"MESSAGE: Writing data to netcdf file: {outfilename}")

        for param in range(0, len(params)):
            if lev == 0:
                # Define netcdf variables from fitting function results
                IO.ncwrite_variables(outfile, [arg_lists[0]["func"].param_names()[param]],
                                     ['f'], ('depth', 'latitude', 'longitude'))

            # Masking function fitting outputs
            params[param].mask = applyMask.create_mask(params[param].mask, [chi_grid],
                                                [-1e10], ['=='], var_look_nan=chi_grid)

            # Adding function fitting parameter to netcdf
            IO.ncwrite_variables(outfile, [arg_lists[0]["func"].param_names()[param]],
                                 [], [], vardata=[params[param]], create_vars=False,
                                 dep_lev=lev)

        # Masking chi_err and obs_err
        obs_err.mask = applyMask.create_mask(obs_err.mask, [chi_grid], [-1e10],
                                             ['=='], var_look_nan=chi_grid)
        chi_grid.mask = applyMask.create_mask(chi_grid.mask, [chi_grid], [-1e10],
                                             ['=='], var_look_nan=chi_grid)

        # Add chi_err and obs_err to netcdf
        IO.ncwrite_variables(outfile, ['obs_err', 'Chi_sq'], [], [],
                             vardata=[obs_err, chi_grid], create_vars=False,
                             dep_lev=lev)

    outfile.close()
