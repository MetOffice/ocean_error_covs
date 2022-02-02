# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
#######################################################################
import numpy as np

class BaseFunction():
    """
        Contains routines that will be common to all function types
    """
    
    def num_params(self):
        """ return number of parameters"""
        return self._num_params
   
 
    def init_guess(self):
        """ return initial guess of parameters"""
        return self._initial_guess

    
    def weights(self):
        """ return weights applied to parameters"""
        return self._weights

    
    def param_names(self):
        """ return names of parameters"""
        return self._param_names

    
    def bounds(self):
        """ return bounds on parameters"""
        return self._bounds

    
    def update_weights(self,new_weights):
        """ Update the values used for the weights"""
        if (len(new_weights) != self._num_params):
            raise ValueError("[ERROR] NEW WEIGHTS IS IN INCORRECT SHAPE")
        self._weights = new_weights

        
    def cost_func(self, func_params, x, y_true):
        """Calculates a cost function
        INPUTS
            func_params:  list of parameters to pass to function
            x:            x values
            y_true:       true values at each x
        """
        return np.sum( (y_true-self.func(x,*func_params))**2. )

        
    def cost_grad(self, func_params, x, y_true):
        """
        Calculates the gradient of the cost functon
         INPUTS
            func_params:  list of parameters to pass to function
            x:            x values
            y_true:       true values at each x       
        """
        grad_func=self.func_grad(x, *func_params)
        grads=[]
        for param, param_grad in zip(func_params, grad_func):
            grads += [np.sum(2*(self.func(x, *func_params) - y_true)*param_grad)]
        return np.array(grads)

    
    def set_lscale_params(self, *args, **kargs):
        """Some derived classes use a lengthscale. 
           To allow for a functional (rather than constant dependence)
           this dummy routine can be overwritten to provide an appropriate 
           interface. Otherwise it does nothing
        """
        return None
    
    
class GaussFunction(BaseFunction):
    """
        Gaussian function
    """
    
    def __init__(self, initial_guess=None, weights=None):
        """
        Initalise function
        INPUTS
            initial_guess:    initial guess for parameters
            weights:         weights to be applies to parameters in minimisation
        """
        self._num_params = 2
        self._initial_guess = initial_guess
        if weights is not None:
            self._weights = weights
        else:
            self.weights = [1.,1.]
        self._param_names = ["Magnitude","Scale"]
        self._bounds = ((0.,None),(0.,None))
        
        
    def func(self, x, mag, lscale):
        """
        Gaussian function
        INPUTS:
            x:       x position for return value
            mag:     magnitude of Guassian
            lscale:  lengthscale of Guassian
        """
        y = mag * np.exp(-(x*x)/(2*lscale*lscale))
        return y
    

    def func_grad(self, x, mag, lscale):
        """
        Gaussian function gradient
        INPUTS:
            x:       x position for return value
            mag:     magnitude of Guassian
            lscale:  lengthscale of Guassian
        """
        xx = x * x
        exp_fctr = np.exp(-(xx)/(2*lscale*lscale))
        dy_dmag = exp_fctr
        dy_dlscale = mag * (xx/lscale**3.) * exp_fctr
        return np.array([dy_dmag,dy_dlscale])
    

    def denormalise_fctr(self, mag_norm, lscale_norm):
        """
        function to return the denormalisation factors to the parameters after 
        fitting
        INPUTS:
            mag_norm:         Normaliation factor used for the magnitude
            lscale_norm:      Normaliation factor used for the length scale
            
        """
        return np.array([mag_norm,lscale_norm])


class Gauss2Function(BaseFunction):
    """
        Sum of 2 guassian functions (optimesd for fitting).
    """
    
    def __init__(self, initial_guess=None, weights=None):
        """
        Initalise function
        INPUTS
            initial_guess:    initial guess for parameters
            weights:         weights to be applies to parameters in minimisation
        """
        self._num_params = 4
        self._initial_guess = initial_guess
        if weights is not None:
            self._weights = weights
        else:
            self.weights = [1.,1.,1.,1.]
        self._param_names = ["Magnitude","Magnitude Ratio",
                             "Long Scale","Scale ratio"]
        self._bounds = ((0.,None),(0.,1.),(0.,None),(0.001,1))
        
        
    def func(self, x, mag, mag_ratio, lscale, lscale_ratio):
        """
        Gaussian2 function
        INPUTS:
            x:               x position for return value
            mag:             magnitude of Function
            mag_ratio:       Ratio in magnitude between the two Guassians
            lscale:          long lengthscales
            lscale_ratio     ratio of the long length scale to the short scale
        """
        y = mag*( mag_ratio*np.exp(-(x*x)/(2*lscale*lscale)) 
                +(1.-mag_ratio)*np.exp(-(x*x)/(2*(lscale_ratio*lscale)**2.)))
        return y
   
 
    def func_grad(self, x, mag, mag_ratio, lscale, lscale_ratio):
        """
        Gaussian function gradient
        INPUTS:
            x:               x position for return value
            mag:             magnitude of Function
            mag_ratio:       Ratio in magnitude between the two Guassians
            lscale:          long lengthscales
            lscale_ratio     ratio of the long length scale to the short scale
        """
        xx = x*x
        exp_fctr1 = np.exp(-(xx)/(2*lscale*lscale))
        exp_fctr2 = np.exp(-(xx)/(2*(lscale_ratio*lscale)**2.))
        dy_dmag = mag_ratio*exp_fctr1 + (1.-mag_ratio)*exp_fctr2
        dy_dmag_ratio = mag * (exp_fctr1-exp_fctr2)
        dy_dlscale = mag*( mag_ratio * (xx/lscale**3.) * exp_fctr1 
                         + (1.-mag_ratio) 
                         * (xx/(lscale_ratio**2.*lscale**3.))*exp_fctr2)
        dy_dlscale_ratio = mag*((1.-mag_ratio) 
                               *(xx/(lscale_ratio**3.*lscale**2.))*exp_fctr2)
        return np.array([dy_dmag, dy_dmag_ratio, dy_dlscale, dy_dlscale_ratio])


    def denormalise_fctr(self, mag_norm, mag_ratio_norm,
                          lscale_norm, lscale_ratio_norm):
        """
        function to return the denormalisation factors to the parameters after 
        fitting
        INPUTS:
            mag_norm:         Normaliation factor used for the magnitude
            mag_ratio_norm:   Normaliation factor used for the ratio
            lscale_norm:      Normaliation factor used for the long length scale
            lscale_ratio:     Ratio of the short length scale to the long            
        """
        return np.array([mag_norm, mag_ratio_norm, lscale_norm, lscale_ratio_norm])
    
    
class MultiGaussFunction(BaseFunction):
        """
            Multi Gaussian function
        """
        
        def __init__(self, initial_guess=None, weights=None, num_funcs=2):
            """
            Initalise function
            INPUTS
                initial_guess:    initial guess for parameters
                weights:         weights to be applied to parameters in 
                                 minimisation
                num_funcs:       number of guassians in function
            """
            self._num_funcs = num_funcs
            self._num_params = 2 * self._num_funcs
            if initial_guess is not None: 
                if len(initial_guess) != self._num_params:
                    raise ValueError("[ERROR] INITIAL GUESS MUST BE NONE OR ALL "
                                     +str(self._num_funcs)+" PARAMETERS MUST BE SET")
            self._initial_guess = initial_guess
                
            if weights is not None:
                if len(weights) == 2*self._num_funcs:
                    self._weights = weights
                else:
                    raise ValueError("[ERROR] WEIGTHS MUST BE SAME LENGTH AS 2*num_func OR NONE")
            else:
                self._weights=[]
                for n in range(self._num_funcs):
                    self._weights += [1.]
                
            self._param_names = []
            self._bounds = []
            for n in range(1, self._num_funcs+1):
                self._param_names += ["Magnitude"+str(n),"Scale"+str(n)]
                self._bounds += [(0.,None),(0.001,None)]

        
        def func(self, x, *params):
            """
            Gaussian function
            INPUTS:
                x:          x position for return value
                params:     array of parameters of the multiGuassian in order 
                            [mag1,lscale1,mag2,lscale2,mag3,lscale3,...]
            """
            #checks
            if len(params) != self._num_params:
                raise ValueError("[ERROR] INCORRECT NUMBER OF MAGNITUDES OR LENGTHSCALES")
            
            y=0
            xx=x*x
            for n in range(0, self._num_params, 2):
                y += params[n] * np.exp(-(xx)/(2.*params[n+1]**2.))
            return y
       
 
        def func_grad(self, x, *params):
            """
            Gaussian function gradient
            INPUTS:
                x:       x position for return value
                mags:       array of magnitudes of the multiGuassian
                params:     array of parameters of the multiGuassian in order 
                            [mag1,lscale1,mag2,lscale2,mag3,lscale3,...]
            """
            #checks
            if len(params) != self._num_params:
                raise ValueError("[ERROR] INCORRECT NUMBER OF MAGNITUDES OR LENGTHSCALES")
            
            xx = x * x
            dy_dparam=[]
            for n in range(0, self._num_params, 2):
            
                exp_fctr = np.exp(-(xx)/(2*params[n+1]**2.))
                dy_dparam += [exp_fctr]
                dy_dparam += [params[n] * (xx/params[n+1]**3.) * exp_fctr]
            
            return np.array(dy_dparam)
        

        def denormalise_fctr(self, mag_norm, lscale_norm):
            """
            function to return the denormalisation factors to the parameters 
            after fitting
            INPUTS:
                mag_norm:         Normaliation factor used for the magnitude
                lscale_norm:      Normaliation factor used for the length scale
                
            """
            return np.array([mag_norm,lscale_norm])


class MultiGaussFunction_FixedLenScale(BaseFunction):
        """
            Multi Gaussian function with fixed langthscales
        """
        
        def __init__(self, initial_guess=None, weights=None, num_funcs=2, lenscales=[]):
            """
            Initalise function
            INPUTS
                initial_guess:    initial guess for parameters
                weights:         weights to be applied to parameters in 
                                 minimisation
                num_funcs:       number of guassians in function
                lenscales:       list lengthscales for each guassian
            """
                
            self._num_funcs = num_funcs
            self._num_params = self._num_funcs
            
            if len(lenscales) != self._num_params:
                raise ValueError("[ERROR] "+str(num_funcs)+" LENGTHSCALES MUST BE SPECIFIED")
            self._lenscales = lenscales
            
            #These variables are used to resolve functional dependencies
            #if any, on the lengthscale parameters
            self._x=None
            self._y=None
            self._norm_scales=1.
            
            if initial_guess is not None: 
                if len(initial_guess) != self._num_params:
                    raise ValueError("[ERROR] INITIAL GUESS MUST BE NONE OR ALL "
                                     +str(self._num_funcs)+" PARAMETERS MUST BE SET")
            self._initial_guess = initial_guess
                
            if weights is not None:
                if len(weights) != self._num_params:
                    raise ValueError("[ERROR] WEIGHTS MUST BE NONE OR ALL "
                                     +str(self._num_funcs)+" PARAMETERS MUST BE SET")
            self._weights = weights
                
            self._param_names=[]
            self._bounds=[]
            for n in range(1, self._num_funcs+1):
                self._param_names += ["Magnitude"+str(n)]
                self._bounds += [(0.,None)]
            self._bounds = tuple(self._bounds)


        def set_lscale_params(self, x, y, norm):
            """Update position"""
            if isinstance(x, int) & isinstance(x, int):
                self._x = x
                self._y = y
                self._norm_scales = norm
            else:
                raise ValueError("[ERROR] X,Y MUST BE INTEGER GRID POSITIONS")
        

        def __resolve_lscale__(self):
            """ 
               Resolves any functions for specifying the lengthscales that 
               depend on x and y
           """
            out_scales=[]
            for n in self._lenscales:
                if callable(n):
                    out_scales+=[n(self._x,self._y)/self._norm_scales]
                else:
                    out_scales+=[n/self._norm_scales]
            return np.array(out_scales)

        
        def func(self, x, *mags):
            """
            Gaussian function
            INPUTS:
                x:          x position for return value
                mags:       magnitudes of the multiGuassian as unnamed arguments
            """
            #checks
            if len(mags) != self._num_funcs:
                raise ValueError("[ERROR] INCORRECT NUMBER OF MAGNITUDES SPECIFIED")
            
            lscale=self.__resolve_lscale__()
            y=0
            xx=x*x
            for n in range(0, self._num_funcs):
                y += mags[n] * np.exp(-(xx)/(2*lscale[n]*lscale[n]))
            return y
        

        def func_grad(self, x, *mags):
            """
            Gaussian function gradient
            INPUTS:
                x:       x position for return value
                mags:       array of magnitudes of the multiGuassian
            """
            #checks
            if len(mags) != self._num_funcs:
                raise ValueError("[ERROR] INCORRECT NUMBER OF MAGNITUDES OR LENGTHSCALES SPECIFIED")
            
            lscale=self.__resolve_lscale__()
            
            xx = x*x
            dy_dmags=[]
            for n in range(0, self._num_funcs):
                dy_dmags += [ 
                    np.exp(-(xx)/(2*lscale[n]*lscale[n])) ]
            return np.array(dy_dmags)
        

        def denormalise_fctr(self, mag_norm):
            """
            function to return the denormalisation factors to the parameters 
            after fitting
            INPUTS:
                mag_norm:         Normaliation factor used for the magnitude
                
            """
            return np.array([mag_norm])
        
        def lengthscales(self):
            """ lengthscales"""
            return self._lenscales
