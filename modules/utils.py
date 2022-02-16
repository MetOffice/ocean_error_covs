# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
from os import path
import numpy as np

class Utils():
      """Class of utility functions"""

      def __init__(self):
          self.fbrmdi = 99999.


      def divide_files_per_proc(self, nproc, files):
          """ Function to divide list of files for each processor

          ***** PARAMETERS ******
          1. nproc: number of processors
          2. files: list of files

          ****** RETURNS ********
          1. list_files: updated list after checking path of all files
          """
          file_lists = [[] for x in range(0, nproc)]
          num_proc=0
          for ffile in files:
              num_proc = num_proc if num_proc < nproc else 0
              file_lists[num_proc] += [ffile]
              num_proc += 1
          return file_lists

      def check_files(self, list_files):
          """ Function to check if files exist

              **** PARAMETERS ****
              1. list_files: list of files

              ***** RETURNS ******
              1. list_files: updated list after checking 
                             path of all files
          """
          checked_files = []
          for archive in list_files:
              if path.exists(archive):
                 checked_files.append(archive)
          if checked_files:
             return checked_files
          else:
             return None


      def dstance(self, lon1, lat1, lon2, lat2):
          """ Returns the distance (in m) between two geographical points.
              Copied from NEMO: NEMO/OPA_SRC/ASM/asminc.F90

          ****** PARAMETERS ******
          1. lon1, lat1: longitude/latitude of geographic point
          2. lon2, lat2: longitude/latitude of another geographic point
                         (or points), needs to be an array

          ******* RETURNS ********
          1. dist: distance between the points in meters (array of same length
                   as lon2 and lat2)
          """
          dls  = np.pi / 180.
          dly1 = lat1 * dls
          dlx1 = lon1 * dls
          dlx2 = lon2[:] * dls
          dly2 = lat2[:] * dls
          dlx = np.sin(dly1) * np.sin(dly2[:]) + np.cos(dly1) * np.cos(dly2[:]) * np.cos(dlx2[:]-dlx1)
          dlx[ np.abs(dlx) > 1. ] = 1.
          dld = np.arctan(np.sqrt(( 1. - dlx )/( 1. + dlx ))) * 222.24 / dls
          dist = dld[:] * 1000.
          return dist


      def map_lon_discontinuity(self, lon, j):
          """ Map the longitude discontinuity for each pair of latitude

          ****** PARAMETERS ******
          1. lon: array with longitudes
          2. j: latitude index

          ******* RETURNS ********
          1. idx1, idx2: i point where the longitude discontinuity
                         occurs for j-1 and j
          """
          diff_lon = np.diff(lon[j-1,:])
          idx1 = np.where(diff_lon == np.amin(diff_lon))[0][0] + 1
          diff_lon = np.diff(lon[j,:])
          idx2 = np.where(diff_lon == np.amin(diff_lon))[0][0] + 1
          return idx1, idx2


      def get_grid_corners(self, lon, lat, i, j, idx, lon_change=[180, -180]):
          """ Get the grid corners

          ****** PARAMETERS ******
          1. lon: array with longitudes
          2. lat: array with latitudes
          3. i, j: i, j indexes
          4. lon_change: the boundaries of the longitude discontinuity

          ******* RETURNS ********
          1. polygon1: containing the grid corners
          2. polygon2: additional corners in case there is a longitude discontinuity
          """
          polygon2 = None
          if i == idx[0] and i == idx[1]:
              polygon1 = np.array([[lon[j-1,i-1], lat[j-1,i-1]], [lon[j,i-1], lat[j,i-1]], \
                                  [lon_change[0], lat[j,i]], [lon_change[0], lat[j-1,i]], \
                                  [lon[j-1,i-1], lat[j-1,i-1]]])
              polygon2 = np.array([[lon_change[1], lat[j-1,i]], [lon_change[1], lat[j,i]], \
                                  [lon[j,i], lat[j,i]], [lon[j-1,i], lat[j-1,i]], \
                                  [lon_change[1], lat[j-1,i]]])
          elif i == idx[0] and idx[0] > idx[1]:
              polygon1 = np.array([[lon_change[1], lat[j-1,i-1]], [lon[j,i-1], lat[j,i-1]], \
                                  [lon[j,i], lat[j,i]], [lon[j-1,i], lat[j-1,i]], \
                                  [lon_change[1], lat[j-1,i-1]]])
          elif i == idx[1] and idx[0] > idx[1]:
              polygon1 = np.array([[lon[j-1,i-1], lat[j-1,i-1]], [lon[j,i-1], lat[j,i-1]], \
                                  [lon_change[0], lat[j,i]], [lon[j-1,i], lat[j-1,i]], \
                                  [lon[j-1,i-1], lat[j-1,i-1]]])
          elif i == idx[0] and idx[0] < idx[1]:
              polygon1 = np.array([[lon[j-1,i-1], lat[j-1,i-1]], [lon_change[1], lat[j,i-1]], \
                                  [lon[j,i], lat[j,i]], [lon[j-1,i], lat[j-1,i]], \
                                  [lon[j-1,i-1], lat[j-1,i-1]]])
          elif i == idx[1] and idx[0] < idx[1]:
              polygon1 = np.array([[lon[j-1,i-1], lat[j-1,i-1]], [lon[j,i-1], lat[j,i-1]], \
                                  [lon[j,i], lat[j,i]], [lon_change[0], lat[j-1,i]],
                                  [lon[j-1,i-1], lat[j-1,i-1]]])
          else:
              polygon1 = np.array([[lon[j-1,i-1], lat[j-1,i-1]], [lon[j,i-1], lat[j,i-1]], \
                                  [lon[j,i], lat[j,i]], [lon[j-1,i], lat[j-1,i]], \
                                  [lon[j-1,i-1], lat[j-1,i-1]]])
          return polygon1, polygon2
