# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:53:31 2015

@author: jonaswallin
"""
import time
import numpy as np

from rtnorm import rtnorm

from rtnorm.oldcode import rtnorm as rntorm_v1

amin = -2.00443204036
amax = 3.48672170399
def test_version1(n = 50000, silent= False):
	"""
		test the speed of rtnorm by Christoph Lassner
		test for both fixed limits a,b and variable limits
		n - number of rep
	"""
	times = np.zeros(3)
	
	start = time.time()
	x = rntorm_v1(amin +0.1, 1,size=n)  # @UnusedVariable
	times[0] = time.time() - start
	
	start = time.time()
	x = rntorm_v1(amax -0.1, 4,size=n)
	times[1] = time.time() - start
	
	
	a_vec = 5*np.random.randn(n)
	b_vec = a_vec + np.abs(5*np.random.randn(n))
	
	start = time.time()
	for i, (a, b) in enumerate(zip(a_vec, b_vec)):
		x[i] =rntorm_v1(a, b)
	times[2] = time.time() - start	
	
	if not silent:
		print('for rntorm_v1:')
		for i, time_ in enumerate(times):
			print('test_{i} = {time:.5f}'.format(i = i, time = time_))
			
	
			
			
	return times
 
def test_version2(n = 50000, silent= False):
	"""
		test the speed of rtnorm through cython
		test for both fixed limits a,b and variable limits
		n - number of rep
	"""
	rtnorm_obj = rtnorm()
	times = np.zeros(4)
	
	start = time.time()
	x = rtnorm_obj.sample(amin +0.1, 1,size=n)  # @UnusedVariable
	times[0] = time.time() - start
	
	start = time.time()
	rtnorm_obj.sample(amax -0.1, 4,size=n)
	times[1] = time.time() - start
	
	
	a_vec = 5*np.random.randn(n)
	b_vec = a_vec + np.abs(5*np.random.randn(n))
	
	start = time.time()
	rtnorm_obj.sample(a_vec, b_vec)
	times[2] = time.time() - start	

	a_vec = 5*np.random.randn(n)
	b_vec = a_vec + np.abs(5*np.random.randn(n))
	mu	 = np.zeros(n)
	sigma  = np.ones(n)
	start = time.time()
	rtnorm_obj.sample(a_vec, b_vec, mu, sigma, error_check=False)
	times[3] = time.time() - start	   
	
	if not silent:
		print('for rntorm_v2:')
		for i, time_ in enumerate(times):
			print('test_{i} = {time:.5f}'.format(i = i, time = time_))
			
	
	
			
			
	return times
	
if __name__ == "__main__":
	
	test_version1()
	test_version2()