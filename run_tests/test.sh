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
scratch="${SCRATCH}/scratch_ocean_error_covs"
#################################################
mkdir -p ${scratch}/HL_error_covs
mkdir -p ${scratch}/HL_function_fitting
mkdir -p ${scratch}/REs
mkdir -p ${scratch}/ancillary

# Set directory and pythonpath
pd=$(pwd)
export PYTHONPATH="${pd}/../:${pd}/run_tests:${PYTHONPATH}"

cp test_files.tar.gz ${scratch}     # copying the files for the test
cp -rf ../KGO ${scratch}            # copying the KGO files for the test
cd ${scratch}

tar -zxf test_files.tar.gz                           # decompress files
mv test_files/model_mesh.nc ${scratch}/ancillary     # moving to ancillary

#################################################################
######## RUNNING SET OF TESTS FOR HL ERROR COVARIANCES ##########
#################################################################
echo $(date)
echo "RUNNING CODE TO CALCULATE HL ERROR COVARIANCES"

cd ${scratch}/HL_error_covs
python3 ${pd}/test_HL_calc.py
if [ $? -ne 0 ]; then
   echo "[ERROR] TEST TO CALCULATE HL ERROR COVARIANCES HAS FAILED"
   echo $(date)
   exit 1
fi

python3 ${pd}/test_valid_results_with_KGOs.py ${scratch}/KGO/HL_error_covs ${scratch}/HL_error_covs "*.nc"
if [ $? -ne 0 ]; then
   echo "[ERROR] KGO TEST FOR HL ERROR COVARIANCES HAS FAILED"
   echo $(date)
   exit 1
fi

#################################################################
cd ${scratch}/HL_function_fitting
cp ${scratch}/HL_error_covs/HL_errorcovs.nc .

echo "RUNNING CODE TO FIT HL ERROR COVARIANCES TO A FUNCTION"
echo $(date)
python3 ${pd}/test_HL_fitting.py
if [ $? -ne 0 ]; then
   echo "[ERROR] TEST TO DO THE FITTING HAS FAILED"
   echo $(date)
   exit 1
fi

python3 ${pd}/test_valid_results_with_KGOs.py ${scratch}/KGO/HL_function_fitting ${scratch}/HL_function_fitting "*.nc"
if [ $? -ne 0 ]; then
   echo "[ERROR] KGO TEST FOR HL ERROR FUNCTION FITTING HAS FAILED"
   echo $(date)
   exit 1
fi

################################################################
####### RUNNING SET OF TESTS FOR REPRESENTATION ERRORS #########
################################################################
echo $(date)
echo "RUNNING CODE TO CALCULATE REPRESENTATION ERRORS"
cd ${scratch}/REs
python3 ${pd}/test_RE.py
if [ $? -ne 0 ]; then
   echo "[ERROR] TEST TO CALCULATE REPRESENTATION ERRORS HAS FAILED"
   echo $(date)
   exit 1
fi

python3 ${pd}/test_valid_results_with_KGOs.py ${scratch}/KGO/REs ${scratch}/REs "*.nc"
if [ $? -ne 0 ]; then
   echo "[ERROR] KGO TEST FOR HL ERROR FUNCTION FITTING HAS FAILED"
   echo $(date)
   exit 1
fi

#################################################################
cd ${pd}
echo "ALL TESTS HAVE BEEN COMPLETED SUCCESSFULLY"
echo $(date)
exit 0
