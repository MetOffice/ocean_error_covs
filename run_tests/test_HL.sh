#!/bin/bash
#################################################
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#################################################
# Enter with scratch folder location to run a simple test
scratch='/net/home/h01/dcarneir/scratch_HL_error_covs'
#################################################
mkdir -p ${scratch}/HL_error_covs
mkdir -p ${scratch}/PostProcessing

# Copying files to scratch folder (error covariance calculation step)
cp ../HL_error_covs/* ${scratch}/HL_error_covs               # copying the code to scratch folder
cp test_HL_calc.py ${scratch}/HL_error_covs                  # copying the script to run the test
cp test_files_HLerrorcov.tar.gz ${scratch}/HL_error_covs     # copying the files for the test

# Copying files to scratch folder (fitting step)
cp ../PostProcessing/* ${scratch}/PostProcessing             # copying the code to scratch folder
cp test_HL_fitting.py ${scratch}/PostProcessing              # copying the script to run the test

######## Calculate HL error covariances with random values for a 2D var ##########
cd ${scratch}/HL_error_covs
tar -zxf test_files_HLerrorcov.tar.gz                        # decompress files 

echo $(date)
echo "RUNNING CODE TO CALCULATE HL ERROR COVARIANCES"
python3 test_HL_calc.py 
if [ $? -ne 0 ]; then
   echo "[ERROR] TEST TO CALCULATE HL ERROR COVARIANCES HAS FAILED"
   echo $(date)
   exit 1
fi

################## Fitting functions to HL error covariances ####################
cd ${scratch}/PostProcessing

# Copying HL results to PostProcessing folder
cp ${scratch}/HL_error_covs/HL_errorcovs.nc .

echo "RUNNING CODE TO FIT HL ERROR COVARIANCES TO A FUNCTION"
echo $(date)
python3 test_HL_fitting.py 
if [ $? -ne 0 ]; then
   echo "[ERROR] TEST TO DO THE FITTING HAS FAILED"
   echo $(date)
   exit 1
fi

echo "TEST HAS BEEN COMPLETED SUCCESSFULLY"
echo $(date)
exit 0
