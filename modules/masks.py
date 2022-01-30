# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np

class applyMask():
      """ Class of masking functions """

      def __init__(self):
         self.undef = 1e10
         self.fbrmdi = 99999.


      def vartype_mask(self, data, source_types, varname='STATION_TYPE'):
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
             return np.ones(len(data.variables[varname][:,0]), dtype=bool)
          
          msk = np.zeros(len(data.variables[varname][:,0]), dtype=bool)
          idobs = data.variables[varname][:,:]
          for id_type in source_types:
              file_type = []
              for i in range(0, idobs.shape[0]):
                  file_type += [int("".join(map(bytes.decode, idobs[i])))]
              file_type = np.array(file_type, dtype=np.int16)
              msk[file_type==id_type] = True

          return msk

      def create_mask(self, mask, vars_ref, target_values, operators,
                      var_look_nan=None, logical='or'):
          """ Function to create mask based on multiple conditions

          ***** PARAMETERS ******
          1. data: data to be masked
          2. vars_ref: list of variables which the masking conditions
                       will be applied
          3. target_values: values which will tell whether masking conditions
                            will be True or False
          4. operators: the logical operators (e.g <, >, or ==)
          5. var_look_nan: variable to look and mask NaNs (othwerise None)
          6. logical: type of logical operation to be performed
                      (e.g. 'or', 'and')

          ********** RETURNS ************
          1. mask: created mask
          """

          for i, target in enumerate(target_values):
              mask = self.mask_operations(mask, vars_ref[i], target,
                             operator=operators[i], logical=logical)

          if var_look_nan is not None:
             mask = np.logical_or(mask, np.isnan(var_look_nan))

          return mask

      def mask_operations(self, basemask, var_ref, target,
                          operator='>', logical='or'):
          """ Function to perform masking operations

          ***** PARAMETERS ******
          1. basemask: mask variable
          2. vars_ref: reference variable which the masking conditions
                       will be applied
          3. target: value which will tell whether masking conditions
                     will be True or False
          4. operator: the logical operator (e.g <, >, ==)
          5. logical: type of logical operation to be performed
                      (e.g. 'or', 'and')

          ********** RETURNS ************
          1. basemask: updated mask
          """

          if logical != 'or' and logical != 'and':
              raise ValueError('[ERROR] MASKING LOGICAL NOT FOUND!')

          if operator == '>':
              if logical == 'or':
                 basemask = np.logical_or(basemask, var_ref > target)
              else:
                 basemask = np.logical_and(basemask, var_ref > target)
          elif operator == '>=':
              if logical == 'or':
                 basemask = np.logical_or(basemask, var_ref >= target)
              else:
                 basemask = np.logical_and(basemask, var_ref >= target)
          elif operator == '<':
              if logical == 'or':
                 basemask = np.logical_or(basemask, var_ref < target)
              else:
                 basemask = np.logical_and(basemask, var_ref < target)
          elif operator == '<=':
              if logical == 'or':
                 basemask = np.logical_or(basemask, var_ref <= target)
              else:
                 basemask = np.logical_and(basemask, var_ref <= target)
          elif operator == '==':
              if logical == 'or':
                 basemask = np.logical_or(basemask, var_ref == target)
              else:
                 basemask = np.logical_and(basemask, var_ref == target)
          elif operator == '!=':
              if logical == 'or':
                 basemask = np.logical_or(basemask, var_ref != target)
              else:
                 basemask = np.logical_and(basemask, var_ref != target)
          else:
              raise ValueError('[ERROR] MASKING CONDITION NOT FOUND!')

          return basemask
