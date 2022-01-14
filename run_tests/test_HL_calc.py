# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
########################################################################################
# Running a simple test case with a 2D random variable to output HL error covariances
import os
import HL_error_covs.master as master
########################################################################################
# Changeable parameters
nproc = 3                                          # number of processors to run the test
filename = 'HL_accumstats.nc'                      # output filename 
bins = [20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800]  # bins of separation distance
inc = 45                                           # how many feedback files for each accum stats output
########################################################################################
# Generate HL accumulated stats files
count = 1
output_files = []
for i in range(1, 90, inc):
    print("MESSAGE: RUNNING SCRIPT TO GENERATE HL ACCUMULATED STATS FILES - STEP {}".format(count))
    list_of_files = ["test_files/HL_random_sample_" + str(n) + ".nc" for n in range(i, i+inc)]
    master.HL_cov_accum_stats(list_of_files, obs_type="VAR", source_types=[34], bins=bins,
                              grid_def=[[-8.0, -2.0, 1.0],[-24.0,-18.0,1.0]],
                              depth_boundaries=[], outfilename=filename, nproc=nproc)
    os.system("mv {} HL_accumstats_{}.nc".format(filename, count))
    output_files.append("HL_accumstats_{}.nc".format(count))
    count = count + 1

# Combine files and calculate HL error covs
master.HL_error_covs(output_files, outfilename='HL_errorcovs.nc')

