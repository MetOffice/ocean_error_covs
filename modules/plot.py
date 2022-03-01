# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import os
import numpy as np
import modules.functions as functions
try:
    if os.environ["HOSTNAME"][0:8] == "expspice":
        import matplotlib
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plotting = True
except:
    plotting = False

class Plots():
      """ Class of plotting functions """

      def __init__(self):
         self.MASK_VAL = -1e10
         self.plotting = plotting

      def plot_data_vs_fitting(self, outfig, plot, x_val, cors, var, obs_err, lats, lons, 
                               depth, params, func_name, num_funcs, lenscale, p_val=None):

          # Creating directory for figures if dont exist already
          os.system("mkdir -p %s " % (outfig))

          # Checking if matplotlib is imported, otherwise skip plot
          if not self.plotting:
             return None
 
          for plot_pos in plot:
              if plot_pos[0]>(len(lons)-1) or plot_pos[1]>(len(lats)-1):
                 continue
              if obs_err[plot_pos[1],plot_pos[0]] == self.MASK_VAL:
                 continue
                
              # Plot lat and lon
              lon = lons[plot_pos[0]]
              lat = lats[plot_pos[1]]
              slat = 'N'
              if(lat<0):
                slat = 'S'
              slon = 'E'
              if(lon<0):
                slon = 'W'

              lscales = []
              if func_name == 'MultiGauss':
                  i_guess = []
                  for n in range(1, len(params), 2):
                      lscales += [params[n][plot_pos[1],plot_pos[0]]]
                      i_guess.append(0.5*var[-1,-1])
                      i_guess.append(params[n][plot_pos[1],plot_pos[0]])
                  plot_func = functions.MultiGaussFunction(initial_guess=i_guess, 
                                              weights=None, num_funcs=num_funcs)

              elif func_name == "MultiGauss_Fixed":
                  for nfunc in range(0, num_funcs):
                      lscales += [lenscale[nfunc][plot_pos[1],plot_pos[0]]] \
                                 if isinstance(lenscale[nfunc], np.ndarray) \
                                 else [lenscale[nfunc]]
                  i_guess = [0.5*var[-1,-1]]*num_funcs
                  plot_func = functions.MultiGaussFunction_FixedLenScale(initial_guess=i_guess, 
                                         weights=None, num_funcs=num_funcs, lenscales=lscales)
              else:
                 raise ValueError('[ERROR] FUNCTION NAME SELECTED IS NOT AVAILABLE')

              # Plot title
              title_string = "Fit %s%s | %s%s | %sm | LengthScales:" % (str(abs(round(lat,2))), 
                                       slat, str(abs(round(lon,2))), slon, str(round(depth,2)))
              for n in range(0, len(lscales)):
                  title_string += " %skm" % (str(round(lscales[n],2)))
              if p_val is not None:
                  title_string += (
                      f" | P value: {p_val[plot_pos[1], plot_pos[0]]:4.3f}")
              f1 = plt.figure()
              plt.title(title_string)

              # Plotting variance
              plt.plot([0], var[plot_pos[1],plot_pos[0]], 'ro', label='Cov. at zero sep. dist.')

              # Plotting covariance
              plt.plot(x_val[:], var[plot_pos[1], plot_pos[0]]*cors[plot_pos[1],plot_pos[0],:], 
                      'bo', label='Covariance')

              # Plotting fitting function
              xx = np.arange(0, x_val[-1]+10., 10.)
              pars=[]
              for param in range(0, len(params)):
                  pars += [params[param][plot_pos[1],plot_pos[0]]]
              plt.plot(xx, plot_func.func(xx, *pars), 'k', label=func_name)
              plt.ylabel('Covariance')
              plt.xlabel('Separation Distance (km)')
              plt.legend()

              plt.savefig("%s/Fit_%s_lat_%s%s_lon_%s%s_dep_%sm.png" % (outfig, func_name,
                          str(abs(round(lat,2))), slat, str(abs(round(lon,2))), slon, str(round(depth,2))))
              plt.close()

