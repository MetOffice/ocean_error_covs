import numpy as np

class applyMask():
      """ Class of masking functions """

      def __init__(self):
         self.undef = 1e10
         self.fbrmdi = 99999.


      def vartype_mask(self, data, source_types):
          """ Mask out those innovations different from observation 
              type defined 

              **** PARAMETERS ****
              1. data: object containing feedback file variables
              2. source_types: Tuple of station types of which innovations 
                               are used for covariance calculations
 
              ****** RETURNS *****
              1. msk: Logical mask array (True for innovations from 
                      required station types)
          """
          if not source_types:
             return np.ones(len(data.variables["LATITUDE"][:]), dtype=bool)
          
          msk = np.zeros(len(data.variables["LATITUDE"][:]), dtype=bool)
          idobs = data.variables["STATION_TYPE"][:,:]
          for id_type in source_types:
              file_type = []
              for i in range(0, idobs.shape[0]):
                  file_type += [int("".join(map(bytes.decode, idobs[i])))]
              file_type = np.array(file_type, dtype=np.int16)
              msk[file_type==id_type] = True
          return msk


      def mask_bad_fdbk_data(self, fdbk_var_array):
          """ Mask bad data from feedback variables

              ***** PARAMETERS *****
              1. fdbk_var_array: object containing arrays of feedback variables
              
              ****** RETURNS *******
              1. lmask: logical mask array
          """ 
          lmask = False
          lmask = np.logical_and(fdbk_var_array.obs_vals[:] != self.fbrmdi,
                                 fdbk_var_array.mod_vals[:] != self.fbrmdi)
          lmask = np.logical_and(lmask, fdbk_var_array.obs_qc[:] == 1)
          return lmask


      def mask_obs_in_grid(self, basemask, ilon, ilat, grid_lon, grid_lat, fdbk_var_array):
          """ Produce a mask for the innovation arrays that masks out observations 
              outside the grid box
      
          ********* PARAMETERS *********
          1. basemask: logical mask array 
          2. ilat, ilon: integer index of the grid box location    
          3. grid_lon, grid_lat: array of longitude and latitude of grid
          4. fdbk_var_arrays: FdbkVarArrays object containing the arrays describing 
                              the innovations
      
          ********** RETURNS ************ 
          1. lmask: logical mask array 
          """
          dlat = np.abs(grid_lat[0] - grid_lat[1])
          dlon = np.abs(grid_lon[0] - grid_lon[1])        
          lmask = np.logical_and(basemask, fdbk_var_array.lons[:] > grid_lon[ilon]-dlon/2.)
          lmask = np.logical_and(lmask, fdbk_var_array.lons[:] <= grid_lon[ilon]+dlon/2.)
          lmask = np.logical_and(lmask, fdbk_var_array.lats[:] > grid_lat[ilat]-dlat/2.)
          lmask = np.logical_and(lmask, fdbk_var_array.lats[:] <= grid_lat[ilat]+dlat/2.)    
          return lmask


      def mask_output_cov_data(self, grid_mean, grid_var, 
            num_obs_in_grid, num_pairs_in_cov, cov_xy, corr_xy):
          """ Mask data prior to output

          ***** PARAMETERS ******
          1. grid mean: Grid mean binned error 
          2. grid var: Grid variance
          3. num_obs_in_grid: Number of obs within squared grid
          4. num_pairs_in_cov: Number of covariance pairs for each bin at each grid
          5. cov_xy: Covariance for each bin of separation distance and at each grid
          6. corr_xy: Correlation for each bin of separation distance and at each grid 
          """
          # account for precision errors by forcing a minimum on variance
          grid_var[grid_var < 0.] = 0.

          grid_mean.mask = np.logical_or(np.ma.abs(grid_mean.mask) > self.undef, np.isnan(grid_mean))
          grid_mean.mask = np.logical_or(grid_mean.mask, num_obs_in_grid == 0)
          grid_var.mask = np.logical_or(np.ma.abs(grid_var) > self.undef, np.isnan(grid_var))
          grid_var.mask = np.logical_or(grid_var.mask, num_obs_in_grid == 0)
          cov_xy.mask = np.logical_or(np.ma.abs(cov_xy) > self.undef, np.isnan(cov_xy))
          cov_xy.mask = np.logical_or(cov_xy.mask, num_pairs_in_cov == 0)
          corr_xy.mask = np.logical_or(np.ma.abs(corr_xy) > self.undef, np.isnan(corr_xy))
          corr_xy.mask = np.logical_or(corr_xy.mask, num_pairs_in_cov == 0)
