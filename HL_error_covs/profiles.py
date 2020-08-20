# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np
import arrays

class ObsProfiles():
      """Class of functions to deal with profile observations"""

      def __init__(self):
         self.undef = 1000


      def random_subsample_profiles(self, depths, depth_range, fdbk_var_array):   
          """ Thin profile observations vertically by choosing a single observation 
              at random within a pre-defined depth range.

              **** PARAMETERS ****
              1. depths: Array of depths of the profile observations
              2. depth_range: Tuple defining boundaries over which to sub-sample
                              profile observations
              3. fdbk_var_arrays: object containing the arrays describing the 
                                  innovations 
              
              ***** RETURNS *****
              1. subsampled_fdbk_var_array: object containing the sub-sampled
                                            arrays describing the innovations
                                            with a single obs for each depth level
          """
          
          # Define new tuples
          lat_dep=[] ; lon_dep =[] ; qc_dep =[] ; mod_dep =[] ; obs_dep = []

          # Loop over obs number
          for n in range(0, len(fdbk_var_array.lats)):
              dep_mask = np.logical_and(depths[n,:] >= depth_range[0],
                                        depths[n,:] <  depth_range[1])
              num_in_mask=np.sum(dep_mask)
              if num_in_mask > 0:
                 lat_dep += [fdbk_var_array.lats[n]]
                 lon_dep += [fdbk_var_array.lons[n]]
                 qc_dep += [fdbk_var_array.obs_qc[n]]
                 random_pick = np.random.random_integers(0, num_in_mask-1)       # NOTE: Results will not be
                 obs_dep += [fdbk_var_array.obs_vals[n,dep_mask][random_pick]]   # reproduceable unless a
                 mod_dep += [fdbk_var_array.mod_vals[n,dep_mask][random_pick]]   # random seed is pre-selected
          
          # define a new fdbk array object
          subsampled_fdbk_var_array = arrays.FdbkVarArrays()
          subsampled_fdbk_var_array.lats = np.array(lat_dep)
          subsampled_fdbk_var_array.lons = np.array(lon_dep)
          subsampled_fdbk_var_array.obs_qc = np.array(qc_dep)
          subsampled_fdbk_var_array.obs_vals = np.array(obs_dep)
          subsampled_fdbk_var_array.mod_vals = np.array(mod_dep)

          # additional masking for profiles
          subsampled_fdbk_var_array.obs_qc[np.abs(subsampled_fdbk_var_array.obs_vals)>self.undef] = self.undef
          return subsampled_fdbk_var_array 
