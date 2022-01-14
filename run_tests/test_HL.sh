#!/bin/bash
#################################################
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#
# Must be run from within ocean_error_covs/run_test
#################################################
# Enter with scratch folder location to run a simple test
scratch="${SCRATCH}/scratch_HL_error_covs"
#################################################
mkdir -p ${scratch}/HL_error_covs
mkdir -p ${scratch}/PostProcessing

# Set directory and pythonpath
pd=$(pwd)
export PYTHONPATH="${pd}/../:${pd}/run_tests:${PYTHONPATH}"

# Copying files to scratch folder (error covariance calculation step)
cp test_files_HLerrorcov.tar.gz ${scratch}/HL_error_covs     # copying the files for the test

######## Calculate HL error covariances with random values for a 2D var ##########
cd ${scratch}/HL_error_covs
tar -zxf test_files_HLerrorcov.tar.gz                        # decompress files 

echo $(date)
echo "RUNNING CODE TO CALCULATE HL ERROR COVARIANCES"
python3 ${pd}/test_HL_calc.py
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
python3 ${pd}/test_HL_fitting.py
if [ $? -ne 0 ]; then
   echo "[ERROR] TEST TO DO THE FITTING HAS FAILED"
   echo $(date)
   exit 1
fi

cd ${pd}
echo "TEST HAS BEEN COMPLETED SUCCESSFULLY"
echo $(date)
exit 0
