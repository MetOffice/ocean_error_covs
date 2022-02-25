# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
########################################################################################
# Running a simple test case with a 2D random variable to calculate representation errors
import os
import representation_errors.RE as RE
########################################################################################
# Changeable parameters
inc = 30                                           # how many feedback files for each RE calculation
thinning = 2
min_num_obs = 4
filename = 'RE.nc'
########################################################################################
# Calculate representation errors for each step
count = 1
output_files = []
model_mesh_file = "../ancillary/model_mesh.nc"
for i in range(1, 90, inc):
    print("MESSAGE: RUNNING SCRIPT TO GENERATE REPRESENTATION ERROR - STEP {}".format(count))
    list_of_files = ["../test_files/HL_random_sample_" + str(n) + ".nc" for n in range(i, i+inc)]
    RE.RE_unresolved_scales(list_of_files, model_mesh_file, obs_type="VAR", source_types=[34],
                            thinning=thinning, min_num_obs=min_num_obs, outfilename=filename,
                            lon_discontinuity=False)
    os.system("mv {} RE_{}.nc".format(filename, count))
    output_files.append("RE_{}.nc".format(count))
    count = count + 1

# Average representation errors
RE.calc_RE_season(output_files, outfilename='RE_final.nc')

