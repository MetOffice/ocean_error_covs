# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
########################################################################################
# Running a simple test case with subset model fields in order to calculate NMC error covariances
import os
from glob import glob
import forecast_diff_error_covs.forecast_diff_error_covs as NMC
########################################################################################
# Changeable parameters
num_cross_covs = 25                                # number of cross covariances
variable = "votemper"                              # variable name
########################################################################################
# Generate NMC model differences
print("MESSAGE: RUNNING SCRIPT TO CALCULATE NMC MODEL DIFFERENCES")
NMC.calc_fld_diffs(sorted(glob("../test_files/model_files_NMC/*.nc")), "./", variable,
                   method="NMC", time_recs=[0, 1], level=None)

inputs = sorted(glob("*diffs.nc"))
num = len(inputs)
input1 = inputs[0:int(num/2)]
input2 = inputs[int(num/2):]
############################################################################
# CONFIG 1 - Horizontal covariances | Num cross covs > 0 | level = 0
print("MESSAGE: CONFIG 1 NMC ERROR COVS | Variable {} | Horizontal covs | Num cross covs = {} | level = 0".format(
      variable, num_cross_covs))

# NMC accumulated stats
NMC.calc_stats_covs(input1, variable, level=0, output_dir='./', stats_file='cov_stats_config1_sample1.nc',
                    num_cross_covs=num_cross_covs, vert=False)
NMC.calc_stats_covs(input2, variable, level=0, output_dir='./', stats_file='cov_stats_config1_sample2.nc',
                    num_cross_covs=num_cross_covs, vert=False)

# NMC combined stats
NMC.combine_stats(sorted(glob("cov_stats_config1*.nc")), variable, output_dir='./',
                  out_file='comb_cov_stats_config1.nc', num_cross_covs=num_cross_covs,
                  vert=False)
#############################################################################
# CONFIG 2 - Horizontal covariances | Num cross covs > 0 | level = 15
print("MESSAGE: CONFIG 2 NMC ERROR COVS | Variable {} | Horizontal covs | Num cross covs = {} | level = 15".format(
      variable, num_cross_covs))

# NMC accumulated stats
NMC.calc_stats_covs(input1, variable, level=15, output_dir='./', stats_file='cov_stats_config2_sample1.nc',
                    num_cross_covs=num_cross_covs, vert=False)
NMC.calc_stats_covs(input2, variable, level=15, output_dir='./', stats_file='cov_stats_config2_sample2.nc',
                    num_cross_covs=num_cross_covs, vert=False)

# NMC combined stats
NMC.combine_stats(sorted(glob("cov_stats_config2*.nc")), variable, output_dir='./',
                  out_file='comb_cov_stats_config2.nc', num_cross_covs=num_cross_covs,
                  vert=False)
###############################################################################
# CONFIG 3 - Horizontal covariances | Num cross covs = 0 | level = 0
print("MESSAGE: CONFIG 3 NMC ERROR COVS | Variable {} | Horizontal covs | Num cross covs = 0 | level = 0".format(
      variable))

# NMC accumulated stats
NMC.calc_stats_covs(input1, variable, level=0, output_dir='./', stats_file='cov_stats_config3_sample1.nc',
                    num_cross_covs=0, vert=False)
NMC.calc_stats_covs(input2, variable, level=0, output_dir='./', stats_file='cov_stats_config3_sample2.nc',
                    num_cross_covs=0, vert=False)

# NMC combined stats
NMC.combine_stats(sorted(glob("cov_stats_config3*.nc")), variable, output_dir='./',
                  out_file='comb_cov_stats_config3.nc', num_cross_covs=0,
                  vert=False)
##############################################################################
# CONFIG 4 - Vertical covariances | Num cross covs > 0
print("MESSAGE: CONFIG 4 NMC ERROR COVS | Variable {} | Vertical covs | Num cross covs = {}".format(
      variable, num_cross_covs))

# NMC accumulated stats
NMC.calc_stats_covs(input1, variable, output_dir='./', stats_file='cov_stats_config4_sample1.nc',
                    num_cross_covs=num_cross_covs, vert=True)
NMC.calc_stats_covs(input2, variable, output_dir='./', stats_file='cov_stats_config4_sample2.nc',
                    num_cross_covs=num_cross_covs, vert=True)

# NMC combined stats
NMC.combine_stats(sorted(glob("cov_stats_config4*.nc")), variable, output_dir='./',
                  out_file='comb_cov_stats_config4.nc', num_cross_covs=num_cross_covs,
                  vert=True)
