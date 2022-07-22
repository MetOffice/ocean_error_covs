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

from multiprocessing import Pool
################## Code modules ##############################
from modules.io_data import IO_netCDF4
from modules.plot import Plots
from modules.masks import applyMask
from modules.fitting import FittingHL
from modules.StatisticTests import StatsTests

# Initialising the classes
IO = IO_netCDF4()
Fitting = FittingHL()
Plots = Plots()
applyMask = applyMask()
StatsTests = StatsTests()


def fitting_function(infile, outfilename, func_name="MultiGauss",
                       num_funcs=2, lenscale=(400,40), plot=None, outfig='./figures', 
                       nproc=4, min_num_obs=2, max_iter=100, scalefac=1.0, f_test=True):

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
    12. f_test perform an f_test on the result and write out the p value
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
    IO.ncwrite_variables(outfile, ['RSS', 'RSS_vs_mean', 'degrees_of_freedom', 'obs_err'],
                         ['f', 'f', 'i', 'f'], ('depth', 'latitude', 'longitude'))
    if f_test:
       IO.ncwrite_variables(outfile, ['P_val'], ['f'], ('depth', 'latitude', 'longitude'))

    # Calculate x positions based on the separation distances
    x_val = Fitting.calc_x_positions(bins)

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
        arg_lists = Fitting.create_arg_list(x_val, cors, var, numobsvar, min_num_obs,
                                            func_name, num_funcs, lenscale, max_iter,
                                            scalefac)
        
        # Get workers to do parallel calculations
        results = workers.map(Fitting.fitter, arg_lists)
        workers.close()

        # Unravel results into output grids
        params, obs_err, rss_func_grid, rss_mean_grid, dof = \
                            Fitting.results_to_grid(results, len(lats), len(lons))

        # If requested perform F-test comparing to mean
        p_val = None
        if f_test:
            if (func_name == "MultiGauss"):
                num_params = 2 * num_funcs
            elif (func_name == "MultiGauss_Fixed"):
                num_params = num_funcs
            else:
                raise ValueError(f"[ERROR] Cannot calculate num_params for function={func_name}")
            p_val = StatsTests.f_test_pvalue(rss_func_grid, rss_mean_grid,
                                             num_params, num_params + dof)

        # Remove scale factor if needed
        obs_err = obs_err/scalefac
        for param in range(0, len(params)):
            if "Magnitude" in arg_lists[0]["func"].param_names()[param]:
               params[param] = params[param]/scalefac

        # Plot some results if requested
        if plot != None:
           print(f"MESSAGE: Plotting results - data versus fitting: {depth[lev]} m")
           Plots.plot_data_vs_fitting(outfig, plot, x_val, cors, var, obs_err, lats, lons,
                                      depth[lev], params, func_name, num_funcs, lenscale, p_val)

        for param in range(0, len(params)):
            if lev == 0:
                # Define netcdf variables from fitting function results
                IO.ncwrite_variables(outfile, [arg_lists[0]["func"].param_names()[param]],
                                     ['f'], ('depth', 'latitude', 'longitude'))

            # Masking function fitting outputs
            params[param].mask = applyMask.create_mask(params[param].mask, [rss_func_grid],
                                                [-1e10], ['=='], var_look_nan=rss_func_grid)

            # Adding function fitting parameter to netcdf
            IO.ncwrite_variables(outfile, [arg_lists[0]["func"].param_names()[param]],
                                 [], [], vardata=[params[param]], create_vars=False,
                                 dep_lev=lev)

        # Masking RSS and obs_err
        rss_func_grid.mask = applyMask.create_mask(rss_func_grid.mask, [rss_func_grid], [-1e10],
                                             ['=='], var_look_nan=rss_func_grid)
        rss_mean_grid.mask = applyMask.create_mask(rss_mean_grid.mask, [rss_func_grid], [-1e10],
                                             ['=='], var_look_nan=rss_func_grid)
        obs_err.mask = applyMask.create_mask(obs_err.mask, [rss_func_grid], [-1e10],
                                             ['=='], var_look_nan=rss_func_grid)

        # Add chi_err and obs_err to netcdf
        IO.ncwrite_variables(outfile, ['obs_err', 'RSS', 'RSS_vs_mean', 'degrees_of_freedom'],
                             [], [], vardata=[obs_err, rss_func_grid, rss_mean_grid, dof],
                             create_vars=False, dep_lev=lev)

        # Add P-val to netcdf
        if f_test:
            p_val.mask = applyMask.create_mask(p_val.mask, [rss_func_grid], [-1e10],
                                              ['=='], var_look_nan=rss_func_grid)
            IO.ncwrite_variables(outfile, ['P_val'], [], [], vardata=[p_val],
                                 create_vars=False, dep_lev=lev)

    outfile.close()
