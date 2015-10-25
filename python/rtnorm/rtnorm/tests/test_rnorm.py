'''
Created on Oct 24, 2015

@author: jonaswallin
'''
import unittest
from rtnorm import rtnorm
import numpy as np
from scipy import stats
#WARNING SCIPY stats.trunc.cdf is numerically unstable!!!!!!!!

#TODO: test mean and variance, density inside the interval
#TODO: test outside side the interval



def KS_test_Rtrunc(rtnorm_obj, a, b, n = 1000, p_lim = 0.01):
	"""
	    a      - limit below
	    b      - limit above
	    n      - number of samples
	    p_lims - limit on which to reject KS 
	"""
	
	f = lambda x: stats.truncnorm.cdf(x, a, b)
	p_val = stats.kstest(rtnorm_obj.sample(a, b,size=n), f)[1]
	if p_val < p_lim:
		print('p_val = {pval}'.format(pval = p_val))
		return False

	return True 



class Test(unittest.TestCase):


	def setUp(self):
		
		self.rtnorm_obj = rtnorm()


	def tearDown(self):
		pass


	def test_inits(self):
		
		self.rtnorm_obj.sample(1)
		self.rtnorm_obj.sample(a = 1)
		self.rtnorm_obj.sample(b = 1)
		self.rtnorm_obj.sample(b = [1., 2.])
		self.assertRaises(ValueError,lambda: self.rtnorm_obj.sample(a = [0.,1.], b = [1.,2.,3.]))
		
	def test_ab(self):
		
		
		self.rtnorm_obj.sample(a = [0.,1.,1.], b = [1.,2.,3.], mu  =1.)
		self.assertRaises(ValueError,lambda: self.rtnorm_obj.sample(a = [0.,2.,1.], b = [1.,2.,3.]))
		
	def test_mu_sigma(self):
		
		self.rtnorm_obj.sample( mu  =np.array([1.,2]))
		self.rtnorm_obj.sample( sigma  =np.array([1.,1]))
		self.assertRaises(ValueError,lambda: self.rtnorm_obj.sample(sigma = [0.,2.,1.],  mu = [1.,2.]))
	
	
	def test_unconstrained(self):
		"""
			if the limits are -inf, inf does it generate samples from  
		"""
		p_val = stats.kstest(self.rtnorm_obj.sample(-np.inf,np.inf,size=1000), 'norm')[1]
		self.assertTrue(p_val > 0.01, msg='KS failed')
		
	def test_random_ab(self):
		n = 10
		a_vec = 10*np.random.randn(n)
		b_vec = a_vec + np.random.rand(n)*10
		for a,b in zip(a_vec, b_vec):
			self.assertTrue(KS_test_Rtrunc(self.rtnorm_obj, a, b),
							msg='KS failed on truncated [{a}, {b}]'.format(a=a, b=b))
		

if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()