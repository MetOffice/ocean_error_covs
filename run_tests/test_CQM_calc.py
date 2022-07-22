# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
########################################################################################
# Running a simple test case with subset model fields in order to calculate CQM error covariances
import os
from glob import glob
import forecast_diff_error_covs.forecast_diff_error_covs as CQM
########################################################################################
# Changeable parameters
num_cross_covs = 25                                # number of cross covariances
variable = "votemper"                              # variable name
########################################################################################
# Generate CQM model differences
print("MESSAGE: RUNNING SCRIPT TO CALCULATE CQM MODEL DIFFERENCES")
CQM.calc_fld_diffs(sorted(glob("../test_files/model_files_CQM/*.nc")), "./", variable,
                   method="CQM", time_recs=[0, 0], level=0)

inputs = sorted(glob("*diffs.nc"))
num = len(inputs)
input1 = inputs[0:int(num/2)]
input2 = inputs[int(num/2):]
############################################################################
# CONFIG 1 - Horizontal covariances | Num cross covs > 0 | level = 0
print("MESSAGE: CONFIG CQM ERROR COVS | Variable {} | Horizontal covs | Num cross covs = {} | level = 0".format(
      variable, num_cross_covs))

# CQM accumulated stats
CQM.calc_stats_covs(input1, variable, level=0, output_dir='./', stats_file='cov_stats_config1_sample1.nc',
                    num_cross_covs=num_cross_covs, vert=False)
CQM.calc_stats_covs(input2, variable, level=0, output_dir='./', stats_file='cov_stats_config1_sample2.nc',
                    num_cross_covs=num_cross_covs, vert=False)

CQM.combine_stats(sorted(glob("cov_stats_config1*.nc")), variable, output_dir='./',
                  out_file='comb_cov_stats_config1.nc', num_cross_covs=num_cross_covs,
                  vert=False)
