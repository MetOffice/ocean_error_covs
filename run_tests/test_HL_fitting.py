# Running a simple test case to fit 2D error covariances to MultiGaussian functions
import master
###################################################################################
# Changeable parameters
nproc = 3                              # Number of processors to run the test
lenscale = (400, 40)                   # Length-scales (MultiGaussian_fixed) or initial guesses (MultiGaussian)
num_funcs = 2                          # Number of Gaussian functions, must match number of lenscales above!
max_iter = 1000                        # Maximum number of iterations
plot_locs = [(1,4), (5,4), (4,1)]      # Location of the plots 
###################################################################################
# Fit HL error covariances to MultiGaussian Function with Fixed Length Scales  
master.HL_fitting_function("HL_errorcovs.nc", "MultiGaussian_fixed_length_scales.nc", 
                           func_name="MultiGauss_Fixed", lenscale=lenscale, plot=plot_locs, 
                           outfig="./figures", nproc=nproc, num_funcs=num_funcs,
                           min_num_obs=20, max_iter=max_iter)

# Fit HL error covariances to MultiGaussian Function
master.HL_fitting_function("HL_errorcovs.nc", "MultiGaussian.nc", 
                           func_name="MultiGauss", lenscale=lenscale, plot=plot_locs, 
                           outfig="./figures", nproc=nproc, num_funcs=num_funcs,
                           min_num_obs=20, max_iter=max_iter)
