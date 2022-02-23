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
import numpy.ma as ma
from multiprocessing import Pool
################## Code modules ##############################
from PostProcessing.io_data import IO
from PostProcessing.posproc import Posproc
from PostProcessing.plot import Plots
from scipy.stats import f

# Initialising the classes
IO = IO()
Posproc = Posproc()
Plots = Plots()

def HL_fitting_function(infile, outfilename, func_name="MultiGauss", 
                       num_funcs=2, lenscale=(400,40), plot=None, outfig='./figures', 
                       nproc=4, min_num_obs=2, max_iter=100, f_test=True):

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
    11. f_test perform an f_test on the result and write out the p value
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
    lats, lons, depth, bins = IO.ncread_dimension_variables(infile)

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename, len(lats), len(lons), len(depth))

    # Write dimension variables (depth, lat, lon and bins)
    IO.ncwrite_dimension_variables(outfile, lats, lons, depth)

    # Add attributes
    outfile.Function = "Function fitting done using the "+func_name+" function"
        
    # Add variables
    IO.nc_define_vars(outfile, "RSS", 'f',
                      ("depth","latitude","longitude"))
    IO.nc_define_vars(outfile, "RSS_vs_mean", 'f',
                      ("depth", "latitude", "longitude"))
    IO.nc_define_vars(outfile, "degrees_of_freedom", 'i',
                      ("depth", "latitude", "longitude"))
    IO.nc_define_vars(outfile, "obs_err", 'f', ("depth","latitude","longitude"))
    if f_test:
        IO.nc_define_vars(outfile, "P_val", 'f',
                          ("depth", "latitude", "longitude"))

    # Calculate x positions based on the separation distances
    x_val = Posproc.calc_x_positions(bins)

    for lev in range(0, len(depth)):
        print("MESSAGE: Fitting function " + func_name + " to ErrorCov data: " + str(depth[lev]) + " m")
        
        # set up workers
        workers = Pool(nproc)
        
        # Reading error covariance variables
        var, cors, numobsvar = IO.ncread_errorcovs(infile, lev)
        
        # Creating list with arguments to run in parallel
        arg_lists = Posproc.create_arg_list(x_val, cors, var, numobsvar, min_num_obs, 
                                            func_name, num_funcs, lenscale, max_iter)
        
        # Get workers to do parallel calculations
        results = workers.map(Posproc.fitter, arg_lists)
        workers.close()

        # Unravel results into output grids
        params, obs_err, rss_func_grid, rss_mean_grid, dof = (
                        Posproc.results_to_grid(results, len(lats), len(lons)))

        # Plot some results if requested
        if plot != None:
           print("MESSAGE: Plotting results - data versus fitting: " + str(depth[lev]) + " m")
           Plots.plot_data_vs_fitting(outfig, plot, x_val, cors, var, obs_err, lats, lons,
                                       depth[lev], params, func_name, num_funcs, lenscale)

        # If requested perform F-test comparing to mean
        if f_test:
            if (func_name == "MultiGauss"):
                num_params = 2 * num_funcs
            elif (func_name == "MultiGauss_Fixed"):
                num_params = num_funcs
            else:
                raise ValueError("Cannot calculate num_params for "
                                 + f"function={func_name}")
            p_val = f_test_pvalue(rss_func_grid, rss_mean_grid, num_params,
                           num_params + dof)
        else:
            p_val = None

        print(f"MESSAGE: Writing data to netcdf file: {outfile}")
        if lev == 0:
           for param in range(0, len(params)):
               # Define netcdf variables from fitting function results
               IO.nc_define_vars(outfile, arg_lists[0]["func"].param_names()[param],
                                 'f', ("depth","latitude","longitude"))

        # Add variables to netcdf
        IO.ncwrite_output(outfile, arg_lists[0]["func"], rss_func_grid,
                          rss_mean_grid, dof, obs_err, params, lev, p_val=p_val)
    
    outfile.close()

def f_test_pvalue(rss_func, rss_mean, num_param, num_points):
    """
    Calculate p_value using an F-test comparing again the fit of the mean.
    Inputs:
        rss_func:   Residual sum of squares against the fitted function
        rss_mean:   Residual sum of squares against the mean
        num_param:  number of parameters for fit
        num_points: number of points in the fit
    Returns:
        p_value of the F-test
    """

    # Calculate the test statistic
    if hasattr(num_points,'mask'):
        msk = (num_points <= 0) | num_points.mask
    else:
        msk = (num_points <= 0)
    dof_numerator = (num_param - 1) * ma.array(np.ones(num_points.shape,
                                                       dtype=int), mask=msk)
    dof_denominator = ma.array(num_points - num_param, mask=msk)
    test_stat = ma.array((((rss_mean - rss_func )/dof_numerator) /
                         (rss_func/dof_denominator)), mask=msk)

    # Calculate the p value
    p_value = ma.array(np.zeros(dof_numerator.shape), mask=msk, dtype=float)
    for n in range(dof_numerator.shape[0]):
        for p in range(dof_numerator.shape[1]):
            if not msk[n, p]:
                f_dist = f(dfn=dof_numerator[n, p], dfd=dof_denominator[n, p])
                p_value[n, p] = f_dist.cdf(test_stat[n, p])

    return p_value
