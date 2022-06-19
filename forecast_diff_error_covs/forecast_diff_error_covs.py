# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
################## Python packages ###########################
import os
from timeit import default_timer as timer
from datetime import timedelta
import numpy as np
import copy as cp
from glob import glob
################## Code modules ##############################
from modules.io_data import IO_xarray
from modules.datasets import StatsDataset

def calc_fld_diffs(list_of_files, output_dir, field_name, level=None,
                   time_recs=[0, 1], method="NMC"):
    """
    Top-level script to generate model field differences based on the NMC 
    or Canadian Quick Method (CQM) method. This is information required 
    to generate error covariance estimates.
    ***** PARAMETERS *****
        1. list_of_files:  List of files containing forecasts
        2. output_dir:     Directory to write the output files to
        3. field_name:     Name of the field to calculate the difference.
        4. level:          Level to calculate differences for (None for all levels)
        5. time_recs:      time indexes to calculate forecast diffs
        6. method:         NMC (National Meteorological Center) or CQM (Canadian Quick Method)
    """

    start_time = timer()
    os.system("mkdir -p {}".format(output_dir))
    model_io = IO_xarray(field_name=field_name, level=level)
    rec2 = None
    for ifile, input_file in enumerate(list_of_files):
        if not os.path.exists(input_file):
            raise ValueError("Model file not found {}".format(input_file))
        data, depthvar = model_io.read_model_files(input_file)
        diff = None

        # Get the data and the time for 1st index in time_recs
        rec1, tim1 = model_io.get_input_model_data(data, time_recs[0])

        # After the first file, calculate the difference between the 
        # time records in the two files.
        if rec2 is not None:

            # NOTE: For the NMC method the time of the two records must be 
            # the same originating from distinct forecasts with same length
            if method=="NMC" and (tim1 - tim2).total_seconds() != 0.0:
                raise ValueError("The times of the two records ",
                                 "are not the same: {} {}".format(tim1, tim2))

            print("Calculating diffs at time {} from file {} and at time {} from file {}".format(
                  tim1, os.path.basename(list_of_files[ifile]),
                  tim2, os.path.basename(list_of_files[ifile-1])))
            diff = cp.deepcopy(rec1)
            if method == "NMC":
               diff.values = np.subtract(rec2.values, rec1.values)
            elif method == "CQM":
               diff.values = np.subtract(rec1.values, rec2.values)
            else:
               raise ValueError("Method not found!")

        # Get the data and the time for the second record in the file
        # for use on the next iteration of the loop over files
        rec2, tim2 = model_io.get_input_model_data(data, time_recs[1])
        if diff is not None:
            output_file = "{}/{}_{}_diffs.nc".format(output_dir,
                                                     tim1.strftime("%Y%m%dT%H%MZ"),
                                                     field_name)
            model_io.write_forecast_errors(data, diff, time_recs[0], output_file)

    end_time = timer()
    print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)))


def calc_stats_covs(input_files, variable, level, output_dir='./',
                    stats_file='cov_stats.nc', num_cross_covs=0,
                    vert=False, wrap=True):
    """
    Top-level routine that calculates accumulated error covariance statistics 
    between pairs of fields from forecast difference.

    ***** PARAMETERS *****
        1. input_files:    List of input files (results from calc_fld_diffs)
        2. variable:       Variable name
        3. level:          Level to calculate error covariances for (none for
                           all levels when calculating vertical covariances)
        4. output_dir:     Output directory
        5.stats_file:      Filename to write the output statistics to
                           If it exists it will be updated
        6.num_cross_covs:  Number of grid points in each N/S/E/W directions
                           to calculate the covariances. Put zero if you don't
                           want to calculate cross covariances.
        7. vert:           Calculate vertical covariances. 
                           If False, calculate horizontal covariances.
        8. wrap:           True => wrap around (global grid). Otherwise don't.
    """
    start_time = timer()
    os.system("mkdir -p %s " % (output_dir))
    print("Calculating error covariances for {} ".format(variable))
    if not vert:
        print("for level {}".format(level))
        model_io = IO_xarray(field_name=variable, level=level)
    else:
        model_io = IO_xarray(field_name=variable)
    if not os.path.exists(input_files[0]):
       raise ValueError("Model file not found {}".format(input_file))

    # Reading list of files containing forecast differences
    fc_error, depthvar = model_io.read_model_files(input_files, concatenate=True,
                                                   concat_dim=model_io.time_dims)
    grid_shape = model_io.get_shape_cov(fc_error, depthvar)
    stats_dataset = StatsDataset(variable, grid_shape, depthvar, num_cross_covs=num_cross_covs,
                                 wrap=wrap, vert=vert)

    stats_dataset.set_up(fc_error[stats_dataset.time_dim])
    for itime in range(fc_error.dims[stats_dataset.time_dim]):
        stats_dataset.add_field(fc_error[variable].isel(time_counter=itime).squeeze())
    stats_dataset.save_stats("{}/{}".format(output_dir, stats_file), save_full=False, overwrite=True)

    end_time = timer()
    print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)))


def combine_stats(input_files, variable, output_dir='/.', out_file='comb_cov_stats.nc',
                  num_cross_covs=0, vert=False, wrap=True):
    """
    Top-level routine that combines error covariance information 
    and calculate derived statistics needed for the function fitting.
    If only one input file is supplied then the output file contains the full
    error covariance information.

    ***** PARAMETERS *****
        1. input_files:    List of input files (results from calc_stats_covs)
        2. variable:       Variable name
        3. output_dir:     Output directory
        4. stats_file:     Filename to write the output statistics to
                           If it exists it will be updated
        5. num_cross_covs: Number of grid points in each N/S/E/W directions
                           to calculate the covariances. Put zero if you don't
                           want to calculate cross covariances.
        6. vert:           Calculate vertical covariances. 
                           If False, calculate horizontal covariances.
        7. wrap:           True => wrap around (global grid). Otherwise don't.
    """
    os.system("mkdir -p %s " % (output_dir))
    print("Combining stats from files {}".format(input_files))

    # Initialising model class
    model_io = IO_xarray(field_name=variable)
    if not os.path.exists(input_files[0]):
       raise ValueError("Model file not found {}".format(input_file))
    start_time = timer()

    # Initialising dataset class
    data, depthvar = model_io.read_model_files(input_files[0])
    in_stats_ds = StatsDataset(variable, None, depthvar, num_cross_covs=num_cross_covs,
                               wrap=wrap, vert=vert)

    # Loop over files
    for stats_file in input_files:
        in_stats_ds.read_stats(stats_file)
        print("Finished reading the input stats file {}".format(stats_file))
        if stats_file == input_files[0]:
            stats_dataset = in_stats_ds
        else:
            stats_dataset.combine_stats(in_stats_ds)
    stats_dataset.save_stats("{}/{}".format(output_dir, out_file), save_full=True)
    end_time = timer()
    print("Elapsed time: {}".format(timedelta(seconds=end_time-start_time)))
