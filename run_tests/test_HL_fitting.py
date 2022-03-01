# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
###################################################################################
# Running a simple test case to fit 2D error covariances to MultiGaussian functions
import HL_error_covs.FunctionFitting as FunctionFitting
import shutil
import netCDF4 as nc4
import numpy as np
###################################################################################
# Changeable parameters
nproc = 3                              # Number of processors to run the test
lenscale = (400, 40)                   # Length-scales (MultiGaussian_fixed) or initial guesses (MultiGaussian)
num_funcs = 2                          # Number of Gaussian functions, must match number of lenscales above!
max_iter = 1000                        # Maximum number of iterations
plot_locs = [(1,4), (5,4), (4,1), (5,5)]      # Location of the plots
###################################################################################
# Fit HL error covariances to MultiGaussian Function with Fixed Length Scales  
FunctionFitting.HL_fitting_function("HL_errorcovs.nc", "MultiGaussian_fixed_length_scales.nc", 
                           func_name="MultiGauss_Fixed", lenscale=lenscale, plot=plot_locs, 
                           outfig="./figures_MGFixed", nproc=nproc, num_funcs=num_funcs,
                           min_num_obs=20, max_iter=max_iter, f_test=False)

# Fit HL error covariances to MultiGaussian Function
FunctionFitting.HL_fitting_function("HL_errorcovs.nc", "MultiGaussian.nc", 
                           func_name="MultiGauss", lenscale=lenscale, plot=plot_locs, 
                           outfig="./figures_MG", nproc=nproc, num_funcs=num_funcs,
                           min_num_obs=20, max_iter=max_iter, f_test=False)

# Fit HL error covariances to MultiGaussian Function with Fixed Length Scales, but with F test conducted
FunctionFitting.HL_fitting_function("HL_errorcovs.nc", "MultiGaussian_fixed_length_scales_Ftest.nc",
                           func_name="MultiGauss_Fixed", lenscale=lenscale, plot=plot_locs,
                           outfig="./figures_MG_Fixed_Ftest", nproc=nproc, num_funcs=num_funcs,
                           min_num_obs=20, max_iter=max_iter, f_test=True)

# Fit HL error covariances to MultiGaussian Function, but with F test conducted
FunctionFitting.HL_fitting_function("HL_errorcovs.nc", "MultiGaussian_Ftest.nc",
                           func_name="MultiGauss", lenscale=lenscale, plot=plot_locs,
                           outfig="./figures_MG_Ftest", nproc=nproc, num_funcs=num_funcs,
                           min_num_obs=20, max_iter=max_iter, f_test=True)

# Repeat the "MultiGaussian Function with F test", but delibratly degrade the
# input to check that the P values decrease

# Set random seed so we always get the same result
random = np.random.default_rng(555)

shutil.copyfile("HL_errorcovs.nc", "HL_errorcovs_degraded.nc")
with nc4.Dataset("HL_errorcovs_degraded.nc",'r+') as degrade_me:

    # Degrade selected points by adding random noise
    shp = degrade_me.variables["Correlation"][:].shape
    degrade_me.variables["Correlation"][0, -3:, -3:, :] += (
        random.normal(0., 0.3, size=(3, 3, shp[-1])))
    degrade_me.variables["Correlation"][:, :, :, :] = (
        np.where(degrade_me.variables["Correlation"][:, :, :, :] > 1.,
                 np.ones(shp), degrade_me.variables["Correlation"][:, :, :, :]))
    degrade_me.variables["Correlation"][:, :, :, :] = (
        np.where(degrade_me.variables["Correlation"][:, :, :, :] < -1.,
                 -1.* np.ones(shp),
                 degrade_me.variables["Correlation"][:, :, :, :]))

# Fit degraqded HL error covariances to MultiGaussian Function with F test.
FunctionFitting.HL_fitting_function("HL_errorcovs_degraded.nc",
                           "MultiGaussian_degraded.nc", func_name="MultiGauss",
                           lenscale=lenscale, plot=plot_locs,
                           outfig="./figures_MG_degraded", nproc=nproc,
                           num_funcs=num_funcs, min_num_obs=20,
                           max_iter=max_iter, f_test=True)
