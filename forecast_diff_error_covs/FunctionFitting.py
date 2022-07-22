# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
################## Python packages ###########################
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
import numpy as np
from timeit import default_timer as timer
from datetime import timedelta
from multiprocessing import Pool
#################### Code modules ############################
from modules.datasets import StatsDataset
from modules.fitting import FittingForecastDiff

def fitting_function(input_stats_file, variable, num_cross_covs=25, vert=False,
                     func_name='MultiGauss', method="fmin_tnc", num_funcs=2,
                     lenscale=(40, 400), min_valid_points=10, nproc=1, output_dir="./",
                     output_file="out_params.nc", anisotropic=True, level=None,
                     max_iter=1000, scalefac=1.0):
    """
    Fit a specified funtion to the error covariances specified in the input_stats_file.

    ***** PARAMETERS *****
    1.  input_stats_file: Input file to read the covariances from
    2.  variable: variable to do the function fitting
    3.  num_cross_covs: number of cross covariances
    4.  vert: true for vertical covariances
    5.  func_name: name of the function
    6.  method: minimisation method
    7.  num_funcs: number of functions
    8.  lenscale: tuple of lengthscales to be applied
    9.  min_valid_points: number of non-masked points to continue with fitting
    10. nproc: number of processors for parallel computation
    11. output_dir: output directory
    12. output_file: output_filename
    13. anisotropic: true for fitting values in distinct directions (e.g. E-W, S-N),
                     false otherwise
    14. level: target model level (None except when vert is True)
    15. max_iter: maximum number of iterations
    16. scalefac: multiplying factor to help convergence when covariances are too small
    """
    # Opening input file
    stats_dataset = StatsDataset(variable, num_cross_covs=num_cross_covs, vert=vert)
    stats_dataset.read_stats(input_stats_file, read_full=True)

    # Initialising function fitting
    fitting = FittingForecastDiff(method=method, func_name=func_name, anisotropic=anisotropic,
                                  num_funcs=num_funcs, lenscale=lenscale, vert=vert, max_iter=max_iter,
                                  min_valid_points=min_valid_points, level=level, scalefac=scalefac)

    # Select variables for fitting
    covs, xvals, lats, lons = fitting.select_vars(stats_dataset.covs_ds)

    # Averaging covariances for each direction
    start_time = timer()
    xvals_avg, covs_avg = fitting.average_covs(xvals, covs)
    end_time = timer()
    print("Elapsed time for averaging: {}".format(timedelta(seconds=end_time-start_time)), flush=True)
    workers = Pool(nproc)

    # Doing function fitting for each direction
    for l in range(0, len(covs_avg)):
        print('Fitting functions in {} directions'.format(fitting.directions[l]), flush=True)
        start_time = timer()
        arg_lists = fitting.create_arg_list(xvals_avg[l], covs_avg[l])
        results = workers.map(fitting.fitter, arg_lists)
        params, chi_grid = fitting.results_to_grid(results, lats.shape[0] , lons.shape[1])
        for v in range(0, len(params)):
            if "var" in fitting.varnames[v]:
               params[v][:,:] = params[v][:,:]/scalefac

        if l > 0:
           params_out = fitting.create_fitting_ds(params, lats, lons, fitting.varout[l],
                                                  fitting.varnames, to_merge=params_out)
        else:
           params_out = fitting.create_fitting_ds(params, lats ,lons,
                                                  fitting.varout[l], fitting.varnames)
        end_time = timer()
        print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)), flush=True)

    workers.close()

    print("Saving the estimated parameters", flush=True)
    params_out.to_netcdf("{}/{}".format(output_dir, output_file))
    end_time = timer()
    print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)), flush=True)
