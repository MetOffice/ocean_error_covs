# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
######################################################################
import numpy as np

class CovPairArrays(): 
    def __init__(self, dims):
        self.x = np.zeros(dims)
        self.y = np.zeros(dims)
        self.x_sq = np.zeros(dims)
        self.y_sq = np.zeros(dims)
        self.xy = np.zeros(dims)
        self.sep_dist = np.zeros(dims)
        
    def append_cov_pair(self, upd):
        self.x = np.append(self.x, upd.x)
        self.y = np.append(self.y, upd.y)
        self.x_sq = np.append(self.x_sq, upd.x_sq)
        self.y_sq = np.append(self.y_sq, upd.y_sq)
        self.xy = np.append(self.xy, upd.xy)
        self.sep_dist = np.append(self.sep_dist, upd.sep_dist)
        return self

        
class CovSumStats():
    def __init__(self, dims):
        self.sum_x = np.ma.zeros(dims)
        self.sum_y = np.ma.zeros(dims)
        self.sum_x_sq = np.ma.zeros(dims)
        self.sum_y_sq = np.ma.zeros(dims)
        self.sum_xy = np.ma.zeros(dims)
        self.num_pairs_in_cov = np.ma.zeros(dims,dtype=np.int32)
        
    def __add__(self, upd):
        update = CovSumStats(self.sum_x.shape)
        update.sum_x = self.sum_x + upd.sum_x
        update.sum_y = self.sum_y + upd.sum_y
        update.sum_x_sq = self.sum_x_sq + upd.sum_x_sq
        update.sum_y_sq = self.sum_y_sq + upd.sum_y_sq
        update.sum_xy = self.sum_xy + upd.sum_xy
        update.num_pairs_in_cov = self.num_pairs_in_cov + upd.num_pairs_in_cov
        return update


class GridSumStats():
    def __init__(self, dims):
        self.num_obs_in_grid = np.ma.zeros(dims, dtype=np.int32)
        self.grid_sum = np.ma.zeros(dims)
        self.grid_sum_obs_std = np.ma.zeros(dims)
        self.grid_sum_sq = np.ma.zeros(dims)
    
    def __add__(self, upd):
        update = GridSumStats(self.grid_sum.shape)
        update.num_obs_in_grid = self.num_obs_in_grid + upd.num_obs_in_grid
        update.grid_sum_obs_std = self.grid_sum_obs_std + upd.grid_sum_obs_std
        update.grid_sum = self.grid_sum + upd.grid_sum
        update.grid_sum_sq = self.grid_sum_sq + upd.grid_sum_sq
        return update


class REstats():
    def __init__(self, dims):
        self.num_obs_in_grid = np.ma.zeros(dims, dtype=np.int32)
        self.RE = np.ma.zeros(dims)
        self.count = np.ma.zeros(dims)

    def __add__(self, upd):
        update = REstats(self.RE.shape)
        update.num_obs_in_grid = self.num_obs_in_grid + upd.num_obs_in_grid
        update.RE = self.RE + upd.RE
        update.count = self.count + upd.count
        return update


class FdbkVarArrays():
    def __init__(self):
        self.obs_vals = ()
        self.obs_std = ()
        self.mod_vals = ()
        self.obs_qc = ()
        self.lons = ()
        self.lats = ()
    
    def __getitem__(self, i):
        update = FdbkVarArrays()
        update.obs_vals = self.obs_vals[i]
        update.obs_std = self.obs_std[i]
        update.mod_vals = self.mod_vals[i]
        update.obs_qc = self.obs_qc[i]
        update.lons = self.lons[i]
        update.lats = self.lats[i]
        return update
