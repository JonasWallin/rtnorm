# -*- coding: utf-8 -*-
"""
Warning this tests can fail by randomness,
should just be ok most of the time....
Created on Fri Oct 23 15:15:06 2015

@author: jonaswallin
"""
from scipy import stats

import numpy as np

from rtnorm.oldcode import rtnorm as rntorm_v1
import unittest



def KS_test_Rtrunc(a, b, n = 1000, p_lim = 0.01):
    """
        a      - limit below
        b      - limit above
        n      - number of samples
        p_lims - limit on which to reject KS 
    """
    
    f = lambda x: stats.truncnorm.cdf(x, a, b)
    p_val = stats.kstest(rntorm_v1(a, b,size=n), f)[1]
    if p_val < p_lim:
        print('p_val = {pval}'.format(pval = p_val))
        return False
        
    return True 

class RtnormTestCase(unittest.TestCase):
    """Tests for variables generated with truncated normal"""

    def test_unconstrained(self):
        """
            if the limits are -inf, inf does it generate samples from  
        """
        p_val = stats.kstest(rntorm_v1(-np.inf,np.inf,size=1000), 'norm')[1]
        self.assertTrue(p_val > 0.01, msg='KS failed')
         
    def test_01(self):
        
        self.assertTrue(KS_test_Rtrunc(0, 1),msg='KS failed on truncated [0,1]')
        
    def test_random_ab(self):
        n = 10
        a_vec = np.random.randn(n)
        b_vec = a_vec + np.random.rand(n)*10
        for a,b in zip(a_vec, b_vec):
            self.assertTrue(KS_test_Rtrunc(a, b),
                            msg='KS failed on truncated [{a}, {b}]'.format(a=a, b=b))
        
if __name__ == '__main__':
    unittest.main()