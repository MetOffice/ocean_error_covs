# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
########################################################################################
# Validating test results with KGOs
import sys
import os
from glob import glob
from netCDF4 import Dataset
import numpy as np

KGO_dir = sys.argv[1]
results_dir = sys.argv[2]
sufix = sys.argv[3]

rtol = 1e-05
########################################################################################
if len(glob(f"{KGO_dir}/{sufix}")) == 0:
    raise ValueError("[ERROR] KGO files not found!")

for f in glob(f"{KGO_dir}/{sufix}"):
    filename = os.path.basename(f)
    nc_kgo = Dataset(f, 'r')
    f_result = f"{results_dir}/{filename}"
    if not os.path.exists(f_result):
       raise ValueError(f"[ERROR] Test results not found for {filename}!")
    nc_results = Dataset(f_result, 'r')

    if len(nc_kgo.variables) != len(nc_results.variables):
       raise ValueError("[ERROR] KGO and test results do not have same number of variables!")

    # Looping over all dimensions from KGO file
    for dim in nc_kgo.dimensions:
        kgo = nc_kgo.dimensions[dim].size
        result = nc_results.dimensions[dim].size
        if kgo != result:
           raise ValueError(f"[ERROR] Dimensions {dim} not matching KGO for {filename}!")

    # Looping over all variables from KGO file
    for var in nc_kgo.variables:
        kgo = nc_kgo.variables[var][:]
        result = nc_results.variables[var][:]
        if not np.allclose(kgo, result, rtol=rtol):
           raise ValueError(f"[ERROR] Variable {var} not matching KGO for {filename}!")
