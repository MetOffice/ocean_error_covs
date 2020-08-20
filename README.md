# Python Error Covariance Estimation

This is a quick User Guide for the **Python3** code which calculates error covariances using [Hollingsworth and Lonnberg](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1600-0870.1986.tb00460.x)'s methodology (see also the same methodology being applied for sea surface temperatures in [Robert-Jones et al. 2016](https://www.sciencedirect.com/science/article/pii/S0034425715302273)). The code is split into two parts:

* **HL_error_covs** (Calculate HL error covariances)
* **PostProcessing** (Fit HL error covariances to a function)

All major parts of the code have detailed documentation strings that give details of their input variables and behavior. To see the documentation strings either look at the code or use the python help command. There is also a folder named **run_tests** where the user can easily run simple tests containing both parts of the code.

For any questions contact Davi Mignac Carneiro (davi.carneiro@metoffice.gov.uk) or Matthew Martin (matthew.martin@metoffice.gov.uk)

## Error Covariance Calculation
The folder **HL_error_covs** contains the following modules:

* **master.py** (top-level functions to calculate HL error covariances)
* **errorCovs.py** (methods involved in the calculation of HL error covariances)
* **io_data.py** (methods to read and write netcdf files)
* **arrays.py** (methods to dynamically create and access python arrays)
* **profiles.py** (methods to deal with profile observations)
* **masks.py** (methods to perform masking operations)
* **utils.py** (utility methods)

Two top-level functions are required to calculate HL error covariances in the script **master.py**. The function **HL_cov_accum_stats** calculates the accumulated error covariance statistics over a specific time period, whereas the function **HL_error_covs** combines multiple files that can be generated by the former function and ultimately calculates the error covariances.

To calculate the accumulated error covariance statistics, the python code needs to be given a list of feedback files (through variable *list_of_fdbackfiles*) and the observation name (through variable *obs_type*). The feedback files should contain the observations and the model interpolated at observation locations (for more details see file examples in **run_tests/test_files_HLerrorcov.tar.gz**). The statistics are then calculated on a grid specified by variable *grid_def* (by default this grid is global and has 1 degree resolution). Horizontal correlations are calculated by averaging statistics into bins at a set distance from each grid point. By default, correlations are calculated every 50 km out to a distance of 1000 km, but this can be changed by specifying a list of bins to the variable *bins*.

Both surface (default) and profile observations can be processed. To process profiles, specify a list of depths into variable *depth_boundaries* which defines the boundaries of the depth levels on which the statistics will be calculated. For surface observations, variable *depth_boundaries* is empty. 

By default the resulting HL error covariances will be written into the netcdf file *corrs.nc*. This can be changed by setting the variable *outfilename*. This output file contains the variances, covariances, and correlations from the Hollingsworth and Lonnberg calculation. It also contains the mean error and the number of observations that went into the calculation.

In some cases observations from multiple sources can exist in one feedback file. In this case, it is possible to process the statistics from a subset of these sources by specifying a list of observation sources in variable *source_types*.

For speed, the code can be run on multiple processors by specifying the variable *nproc*.

## Post-Processing

The post-processing takes the HL error covariance netcdf file generated by the statistics part of the code (described above) and fits a user-specified function to this data. The folder **PostProcessing** contains the following modules:

* **master.py** (top-level function to do the fitting)
* **posproc.py** (methods to perform all the post-processing steps)
* **functions.py** (methods to invoke the fitting functions)
* **io_data.py** (methods to read and write netcdf files)
* **plot.py** (methods to plot fitting results)

The function to be fitted has to be specified in variable *func_name*. At present the only options are a Multi-Gaussian (option *MultiGauss*) or a Multi-Gaussian with fixed-length scales (option *MultiGauss_Fixed*). These functions are defined in **functions.py** and can have an arbitrary number of Gaussian components (default 2, set with variable *num_funcs*). For the Multi-Guassian with fixed-length scales, the length scales are specified as a tuple given to variable *lenscale*. For the Multi-Gaussian function, the variable *lenscale* is also used as the function's initial guess. The maximum number of iterations for the fitting is defined by the variable *max_iter*. 

Output from the function fitting is written into a user-specified file, containing the chi-squared values, observational errors, magnitude of the background error variances and their length-scales for each grid point. Plots comparing the original covariances and the function best fit can also be generated for specific locations by providing the grid positions of X,Y pairs to variable *plot*. If no plots are needed simply set variable *plot* to False. As with the statistics code, the post-processing code can be run on multiple processors by specifying variable *nproc*.

It is important to note that additional post-processing may also be necessary, such as smoothing the fitted data, which is not currently included in the code.

## Requirements to run
* Python3 installed 
* Package [netCDF4](https://unidata.github.io/netcdf4-python/netCDF4/index.html) for python
* [NumPy](https://pypi.org/project/numpy/) (code validated with version 1.15.4)
* [SciPy](https://www.scipy.org/) (code validated with version 1.2.1)
* [Matplotlib](https://matplotlib.org/) in case plot options are requested by the user

## Running a simple test

Feedback files for a 2D variable with random values were generated for a small subset of the global domain and are stored in **run_tests/test_files_HLerrorcov.tar.gz**. The instructions to run the test are given below: 

* In the folder **run_tests**, the user should enter the path where the test will be performed by changing the variable *scratch* in the script **test_HL.sh**. 

* After setting the path then just run the test with the command **bash test_HL.sh**.

* The shell script will call two python scripts: **test_HL_calc.py** to calculate the HL error covariances and **test_HL_fitting.py** to fit the error covariances to a function.
  
* All the netcdf files generated by the test are stored in the path provided by the user, as well as the plots for specific locations.

* The user can play with changing some input parameters of the test. In each python script there is a section of changeable parameters where the user can change the number of processors, the bins of separation distance, plot locations, etc.

## License

[BSD-3 License](LICENSE)
