# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ocean_error_covs and is released under the BSD 3-Clause
# license.
# See LICENSE in the root of the repository for full licensing details.
################### Python packages ###########################
import modules.functions as function
import numpy as np
import numpy.ma as ma
from scipy.stats import f

class StatsTests():
    """
    Statistical tests that can be applied to calculated fits
    """

    def __init__(self):
        pass

    def f_test_pvalue(self, rss_func, rss_mean, num_param, num_points):
        """
        Calculate p_value using an F-test comparing against the fit of the mean.
        Inputs:
            1. rss_func:   Residual sum of squares against the fitted function
            2. rss_mean:   Residual sum of squares against the mean
            3. num_param:  number of parameters for fit
            4. num_points: number of points in the fit
        Returns:
            1. p_value of the F-test
        """
    
        # Calculate the test statistic
        if hasattr(num_points, 'mask'):
            msk = (num_points <= 0) | num_points.mask
        else:
            msk = (num_points <= 0)
        dof_numerator = (num_param - 1) * ma.array(np.ones(num_points.shape,
                                                           dtype=int), mask=msk)
        dof_denominator = ma.array(num_points - num_param, mask=msk)
        test_stat = ma.array((((rss_mean - rss_func)/dof_numerator) /
                             (rss_func/dof_denominator)), mask=msk)
    
        # Calculate the p value
        p_value = ma.array(np.zeros(dof_numerator.shape), mask=msk, dtype=float)
        for n in range(dof_numerator.shape[0]):
            for p in range(dof_numerator.shape[1]):
                if not msk[n, p]:
                    f_dist = f(dfn=dof_numerator[n, p], dfd=dof_denominator[n, p])
                    p_value[n, p] = f_dist.cdf(test_stat[n, p])
    
        return p_value

