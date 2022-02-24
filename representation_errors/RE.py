# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
################## Python packages ###########################
import os
import numpy as np
################## Code modules ##############################
import modules.arrays as arrays
from modules.io_data import IO
from modules.utils import Utils
from modules.masks import applyMask
from representation_errors.calcREs import CalcREs

# Initialising the classes
IO = IO()
Utils = Utils()
applyMask = applyMask()
CalcREs = CalcREs()

def RE_unresolved_scales(list_of_fdbackfiles, model_mesh_file, obs_type="SLA",
                 outfilename="RE_stats.nc", source_types=[], qc_val=[], min_num_obs=5,
                 depth_lev=[0], latcen='gphit', loncen='glamt', latcor='gphit',
                 loncor='glamf', mask='tmaskutil', thinning=1, lon_discontinuity=True):

    """ Top-level routine that calculates the representation error due to
        model unresolved scales based on Oke & Sakov (2008).
        It calculates the average of all observations falling within
        a grid cell and then computes their standard deviation. These
        stats are then summed and stored.

    *************** PARAMETERS *******************
    1. list_of_fdbackfiles: list of feedback files
    2. model_mesh_file: model mesh filename
    3. obs_type: observation type to process
    4. outfilename: name of the output file
    5. source_types: observation id types to process
                     (default: [], which means process all)
    6. qc_val: QC values to filter obs
                     (default: [], which means do not filter any obs)
    7. min_num_obs: minimum number of obs to proceed with calculations
    8. depth: list with depth levels (default is [0] - only surface)
    9. latcen: variable name of the latitudes at the center of grid cell
    10. loncen: variable name of the longitudes at the center of grid cell
    11. latcor: variable name of the latitudes at the corners of grid cell
    12. loncor: variable name of the longitudes at the corners of grid cell
    13. mask: mask array (None in case mask is not provided)
    14. thinning: number to subsample obs
    15. lon_discontinuity: True in case there is discontinuity in longitude
                           (e.g. 180 to -180)
    """
    # Check list of feedack files
    list_of_fdbackfiles = Utils.check_files(list_of_fdbackfiles)
    if list_of_fdbackfiles is None:
       raise ValueError("[ERROR] FEEDBACK FILES NOT FOUND")

    # Check model mesh file
    model_mesh_file = Utils.check_files([model_mesh_file])
    if list_of_fdbackfiles is None:
       raise ValueError("[ERROR] MODEL_MESH_FILE NOT FOUND")

    latmod = IO.ncread_variables(model_mesh_file[0], [latcen], dep_lev=0)[0]
    lonmod = IO.ncread_variables(model_mesh_file[0], [loncen], dep_lev=0)[0]

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename, ['y', 'x', 'z'],
                                     [latmod.shape[0], lonmod.shape[1], len(depth_lev)])

    # Write dimension variables
    IO.ncwrite_variables(outfile, ['latitude'], ['f'], ('y', 'x'), [latmod])
    IO.ncwrite_variables(outfile, ['longitude'], ['f'], ('y', 'x'), [lonmod])
    IO.ncwrite_variables(outfile, ['depth'], ['f'], ('z'), [depth_lev])

    # Create 3D netcdf variables (accumulated stats)
    IO.ncwrite_variables(outfile, ['RE', 'GridNumObs', 'Count'],
                         ['f', 'i', 'i'], ('z', 'y', 'x'))

    for dep_lev in range(0, len(depth_lev)):
        print("MESSAGE: Calculating accumulated RE stats at {} m".format(depth_lev[dep_lev]))

        arg_list = {"list_of_files": sorted(list_of_fdbackfiles),
                    "model_file": model_mesh_file[0],
                    "obs_type": obs_type,
                    "source_types": source_types,
                    "qc_val": qc_val,
                    "min_num_obs": min_num_obs,
                    "latcor": latcor,
                    "loncor": loncor,
                    "mask": mask,
                    "depth": depth_lev[dep_lev],
                    "thinning": thinning,
                    "lon_discontinuity": lon_discontinuity}

        # Calculating REs for a specific time window
        RE = CalcREs.calc_RE_unresolved_scales(arg_list)

        # Write RE stats to output file
        IO.ncwrite_variables(outfile, ['RE', 'GridNumObs', 'Count'], [], [],
                             vardata=[RE.RE, RE.num_obs_in_grid , RE.count],
                             create_vars=False, dep_lev=dep_lev)
    outfile.close()


def calc_RE_season(list_of_files, outfilename="RE.nc"):

    """ Top-level routine that averages the Representation Errors from a
        list of netcdf files containing the RE for short-time windows

    *************** PARAMETERS *******************
    1. list_of_files: list of netcdf files containing the accumulated statistics
    2. outfilename: name of file to contain the output statistics
    """
    # Check list of netcdf files
    list_of_files = Utils.check_files(list_of_files)
    if list_of_files is None:
       raise ValueError("[ERROR] NETCDF FILES NOT FOUND")

    ncdata = IO.ncread_variables(list_of_files[0], ['latitude', 'longitude', 'depth'])
    lat = ncdata[0]
    lon = ncdata[1]
    depths = ncdata[2]

    # Create netcdf object and add dimensions
    outfile = IO.nc_define_dimensions(outfilename, ['y', 'x', 'z'],
                                      [lat.shape[0], lon.shape[1], depths.shape[0]])

    # Write dimension variables
    IO.ncwrite_variables(outfile, ['latitude'], ['f'], ('y', 'x'), [lat])
    IO.ncwrite_variables(outfile, ['longitude'], ['f'], ('y', 'x'), [lon])
    IO.ncwrite_variables(outfile, ['depth'], ['f'], ('z'), [depths])

    # Add netcdf variables
    IO.ncwrite_variables(outfile, ['RE', 'GridMeanNumObs'], ['f', 'f'],
                         ('z', 'y', 'x'), fill_value=[1e10, 1e10])


    for dep_lev in range(0, depths.shape[0]):
        print("MESSAGE: Calculating RE at {} m".format(depths[dep_lev]))
        final_RE = arrays.REstats((lat.shape[0], lon.shape[1]))
        for f in list_of_files:

            # Read RE variables from netcdf
            ncdata = IO.ncread_variables(f, ['RE', 'GridNumObs', 'Count'], dep_lev=dep_lev)
            RE = arrays.REstats((lat.shape[0], lon.shape[1]))
            RE.RE = ncdata[0]
            RE.num_obs_in_grid = ncdata[1]
            RE.count = ncdata[2]

            final_RE += RE

        # Averaging RE over a season
        RE = np.where(final_RE.count > 0, final_RE.RE/final_RE.count, 0)
        GridMeanNumObs = np.where(final_RE.count > 0, final_RE.num_obs_in_grid/final_RE.count, 0)

        # Masking output data
        msk = applyMask.create_mask(np.zeros((lat.shape[0], lon.shape[1])),
                                   [final_RE.count], [0], ['<='])
        RE[msk] = 1e10
        msk = applyMask.create_mask(np.zeros((lat.shape[0], lon.shape[1])),
                                   [final_RE.count], [0], ['<'])
        GridMeanNumObs[msk] = 1e10

        # Write variables to output file
        IO.ncwrite_variables(outfile, ['RE', 'GridMeanNumObs'], [], [],
                             vardata=[RE, GridMeanNumObs],
                             create_vars=False, dep_lev=dep_lev)

    outfile.close()
