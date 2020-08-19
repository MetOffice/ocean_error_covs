import numpy as np
import numpy.ma as ma
import scipy.optimize as opt
import functions

class Posproc():
      """ Class of post-processing functions """

      def __init__(self):
         self.MASK_VAL = -1e10  


      def calc_x_positions(self, bins):
          """ Calculate x positions based on the bins of separation distance

          ******** PARAMETERS **********
          1. bins: bins of separation distance

          ******** RETURNS *************
          1. x_val: values of x positions
          """
          x_val=[0.5*bins[0]]
          for n in range(0, len(bins)-1):
              x_val += [0.5 * (bins[n] + bins[n+1])]
          x_val = np.array(x_val)
          return x_val


      def create_arg_list(self, x_val, cors, var, numobsvar, min_num_obs, 
                          func_name, num_funcs, lenscale, max_iter):
          """ Create list of arguments for running in parallel

          ******** PARAMETERS **********
          1. xval: list of x positions
          2. cors: correlation of innovation statistics
          3. var: grid box variance
          4. numobsvar: number of obs within grid box
          5. min_num_obs: minimum number of obs to carry on with fitting
          6. func_name: name of fitting function
          7. num_funcs: number of functions to use
          8. lenscale: lengthscales to apply
          9. max_iter: max number of iterations

          ******** RETURNS *************
          1. arg_list: list with arguments
          """
          arg_lists = [] 
          for nn in range(0, cors.shape[0]):
              for pp in range(0, cors.shape[1]):

                  if numobsvar[nn,pp] < min_num_obs:
                      cors[nn,pp,:] = np.nan
                      var[nn,pp] = np.nan

                  arg_lists += [{ "covs": cors[nn,pp,:]*var[nn,pp], # NOTE: Reweighting by the central point
                                  "x_vals": x_val[:],
                                  "var": var[nn,pp],
                                  "max_iter": max_iter }]

                  lscales = []
                  for nfunc in range(0, num_funcs):
                      lscales += [lenscale[nfunc][nn,pp]] if isinstance(lenscale[nfunc], np.ndarray) \
                                                         else [lenscale[nfunc]]

                  if func_name == 'MultiGauss':
                      i_guess = []
                      for i in range(0, num_funcs):
                          i_guess.append(0.5*var[nn,pp])
                          i_guess.append(lscales[i])
                      arg_lists[-1]["func"] = functions.MultiGaussFunction(initial_guess=i_guess, 
                                                              weights=None, num_funcs=num_funcs)
                  elif func_name == 'MultiGauss_Fixed':
                      i_guess = [0.5*var[nn,pp]] * num_funcs
                      arg_lists[-1]["func"] = functions.MultiGaussFunction_FixedLenScale(initial_guess=i_guess, 
                                                        weights=None, num_funcs=num_funcs, lenscales=lscales)
                  else:
                    raise ValueError('[ERROR] FUNCTION NAME SELECTED IS NOT AVAILABLE')

          return arg_lists


      def fitter(self, args):
          """
          Algorithm to do the fitting

          ************ PARAMETERS ***********
          1. args: dictionary containing correlations and variance for 1 point,
                   instance of fitting function and x values

          ************* RETURNS *************
          1. fit: fitting results
          2. obs_err: obs error 
          3. chi: chi squares criterion to check quality of fitting procedure
          """
          
          # Read arguments
          covs = args["covs"]
          func = args["func"]
          x_vals = args["x_vals"]
          var = args["var"]
          max_iter = args["max_iter"]
              
          #check the corrs for nans and calculate fit
          nans = np.isnan(covs) 
          if np.sum(nans) == covs.shape[0]:
              fit = self.MASK_VAL*np.ones((func.num_params()))
              obs_err = self.MASK_VAL
              chi = self.MASK_VAL
          else:
              #NOTE: we do not use the zero point in the fitting in order to get the obs error
              fit = opt.fmin_tnc(func.cost_func, func.init_guess(), disp=0, fprime=func.cost_grad,
                                 maxfun=max_iter, args=(x_vals[:], covs[:]),
                                 bounds=func.bounds())[0]
              obs_err = var - func.func(0, *fit)
              chi = func.cost_func(fit, x_vals[:], covs[:])
          return (fit, obs_err, chi)


      def results_to_grid(self, results, nlat, nlon):
          """ Grid the results from multiprocessing 

          ************ PARAMETERS ***********
          1. results: results merged over all the processors
          2. nlat: number of latitudes
          3. nlon: number of longitudes

          ************* RETURNS *************
          1. params: fitting parameters gridded
          2. obs_err: obs error gridded
          3. chi: chi squares gridded
          """

          chi_grid=ma.zeros((nlat, nlon))
          obs_err=ma.zeros((nlat, nlon))
          num_params = len(results[0][0][:])
          params = []
          for param in range(0, num_params):
              params += [ma.zeros((nlat, nlon))]

          n = 0
          for nn in range(0, nlat):
              for pp in range(0, nlon):
                for param in range(0, num_params):
                    params[param][nn,pp] = results[n][0][param]
                obs_err[nn,pp] = results[n][1]
                chi_grid[nn,pp] = results[n][2]
                n += 1

          return params, obs_err, chi_grid
