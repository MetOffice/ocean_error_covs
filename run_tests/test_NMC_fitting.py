# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
###################################################################################
# Running simple test cases to fit NMC error covariances to various functions
import forecast_diff_error_covs.FunctionFitting as FunctionFitting
###################################################################################
# Changeable parameters
variable = "votemper"                  # Variable name
nproc = 3                              # Number of processors to run the test
max_iter = 1000                        # Maximum number of iterations
###################################################################################
# CONFIG 1 - Horizontal anisotropic covs | MultiGauss | min method = fmin_tnc
print("MESSAGE: CONFIG 1 NMC FUNC FIT | Horizontal anisotropic covs | MultiGauss | min method = fmin_tnc")
FunctionFitting.fitting_function('comb_cov_stats_config1.nc', variable, num_cross_covs=25,
                                 func_name='MultiGauss', method="fmin_tnc", num_funcs=2,
                                 lenscale=(40, 400), nproc=nproc, output_file="fit_config1.nc", 
                                 anisotropic=True, max_iter=max_iter)

# CONFIG 2 - Horizontal anisotropic covs | MultiGauss | min method = fmin_tnc | scalefac = 10
print("MESSAGE: CONFIG 2 NMC FUNC FIT | Horizontal anisotropic covs | MultiGauss | min method = fmin_tnc | scalefac = 10")
FunctionFitting.fitting_function('comb_cov_stats_config1.nc', variable, num_cross_covs=25,
                                 func_name='MultiGauss', method="fmin_tnc", num_funcs=2,
                                 lenscale=(40, 400), nproc=nproc, output_file="fit_config2.nc", 
                                 anisotropic=True, max_iter=max_iter, scalefac=10.0)

# CONFIG 3 - Horizontal isotropic covs | MultiGauss | min method = fmin_tnc
print("MESSAGE: CONFIG 3 NMC FUNC FIT | Horizontal isotropic covs | MultiGauss | min method = fmin_tnc")
FunctionFitting.fitting_function('comb_cov_stats_config1.nc', variable, num_cross_covs=25,
                                 func_name='MultiGauss', method="fmin_tnc", num_funcs=2,
                                 lenscale=(40, 400), nproc=nproc, output_file="fit_config3.nc", 
                                 anisotropic=False, max_iter=max_iter)

# CONFIG 4 - CONFIG NMC FUNC FIT | Horizontal anisotropic covs | MultiGauss | min method = L-BGFS-B"
print("MESSAGE: CONFIG 4 NMC FUNC FIT | Horizontal anisotropic covs | MultiGauss | min method = L-BGFS-B")
FunctionFitting.fitting_function('comb_cov_stats_config1.nc', variable, num_cross_covs=25,
                                 func_name='MultiGauss', method="L-BFGS-B", num_funcs=2,
                                 lenscale=(40, 400), nproc=nproc, output_file="fit_config4.nc", 
                                 anisotropic=True, max_iter=max_iter)

# CONFIG 5 - Horizontal anisotropic covs | MultiGauss Fixed | min method = fmin_tnc
print("MESSAGE: CONFIG 5 NMC FUNC FIT | Horizontal anisotropic covs | MultiGauss Fixed | min method = fmin_tnc")
FunctionFitting.fitting_function('comb_cov_stats_config1.nc', variable, num_cross_covs=25,
                                 func_name='MultiGauss_Fixed', method="fmin_tnc", num_funcs=2,
                                 lenscale=(40, 400), nproc=nproc, output_file="fit_config5.nc", 
                                 anisotropic=True, max_iter=max_iter)

# CONFIG 6 - Vertical anisotropic covs | Gauss | min method = fmin_tnc
print("MESSAGE: CONFIG 6 NMC FUNC FIT | Vertical anisotropic covs | Gauss | min method = fmin_tnc")
FunctionFitting.fitting_function('comb_cov_stats_config4.nc', variable, num_cross_covs=25,
                                 func_name='Gauss', level=0, method="fmin_tnc", num_funcs=1,
                                 lenscale=(10,), nproc=nproc, output_file="fit_config6.nc", 
                                 anisotropic=True, max_iter=max_iter, vert=True)
