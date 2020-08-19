**** Test to calculate HL error covariances and fit them to a function *****
****************************************************************************

1) Feedback files for a 2D variable with random values were generated for a small subset of the domain 
   (test_files_HLerrorcov.tar.gz)

2) In order to run the test you need first to provide the path where the test folder will be located by changing 
   the variable "scratch" in the script "test_HL.sh". After setting the path, just run "bash test_HL.sh".

3) In the first step of the test, the HL error covariances will be calculated by running the python3 script
   "test_HL_calc.py". This is done in two tasks: the first one generates multiple files of accumulated 
   error covariance statistics and the second one combines these files in order to calculate the HL error 
   covariances for the whole period. 

4) The netcdf file containing the HL error covariances are then used by the python3 script "test_HL_fitting.py",
   which will fit the error covariances to a MultiGaussian function. These functions could be either a MultiGaussian
   function or a MultiGaussian function with pre-defined length scales. 

5) All the netcdf files generated are stored in the test folder as well as the figures showing the function fitting 
   for specific locations over the small domain subset. 

6) You can play with changing the test parameters. In each python script there is a section of "changeable parameters"
   where you can change the number of processors, the bins of separation distance, plot locations, etc. 

In case of any questions please contact Davi Carneiro (davi.carneiro@metoffice.gov.uk)

I hope you have good fun! :)
