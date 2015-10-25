'''
Created on Oct 24, 2015

@author: jonaswallin
'''

from .cython.rtnorm import rtnormClass as cyrtnormClass  # @UnresolvedImport
import numpy as np

from scipy.special import erf



class rtnormClass(object):
	'''
		Wrapper class for cython object
		and check validity of inputs
	'''
	def __init__(self):
		
		
		self.cyrtnormClass = cyrtnormClass()
		
	

	def error_checking(self, a, b, mu, sigma, size):
		'''
			Helper fucntion for sample make sures that things are correct with the parameters
		
		'''
		shapes = []
		max_shape = 1
		
		a	 = np.array(a)
		b	 = np.array(b)
		mu	= np.array(mu)
		sigma = np.array(sigma)
		list_obj  = [a,b,mu,sigma]
		for i in range(len(list_obj)):
			if list_obj[i].shape == ():
				list_obj[i] = np.array([list_obj[i]])
			if len(list_obj[i].shape ) > 2:
				raise ValueError('all input must be of dim 1')
			shapes.append(list_obj[i].shape[0])
			
			if max_shape == 1:
				max_shape = list_obj[i].shape[0]
				
			if list_obj[i].shape[0] > 1 and list_obj[i].shape[0] != max_shape and max_shape != 1:
				raise ValueError('all input must be of equal dimension or dimension one')
		
		
		
		if max_shape > 1:
			
			if size is not None:
				raise ValueError('size can only be used when all input is scalars')
			
			for i in range(len(list_obj)):
				if list_obj[i].shape[0] != max_shape:
					list_obj[i] = np.repeat(list_obj[i], max_shape)


		else:
			
			if size is not None:
				
				for i in range(len(list_obj)):
					list_obj[i] = np.repeat(list_obj[i], size)
					
		for i in range(len(list_obj)):
			list_obj[i] = list_obj[i].astype(np.double)
		
		
		if np.all(list_obj[0] < list_obj[1]) == False:
			raise ValueError('a should be smaller then b')
		
		return list_obj[0], list_obj[1], list_obj[2], list_obj[3]
		
	def sample(self, a		   = np.array([-np.Inf])
				   , b		   = np.array([np.Inf])
				   , mu		  = np.array([0.])
				   , sigma	   =  np.array([1.])
				   , size		= None
				   , error_check = True):
		"""
			Samples from truncated normal in the interval [a,b]
			
			*a*		   - (nx1) left limits
			*b*		   - (nx1) right limits
			*mu*		  - (nx1) expectation
			*sigma*	   - (nx1) the standard devations
			*size*		  - (int) number of samples only applicaple when things all params are in non vector form
			*error_check* - do you know what you are doing?
		"""
		
		if error_check:
			a, b, mu, sigma = self.error_checking(a, b, mu, sigma, size)
			
		
		X = self.cyrtnormClass.sample(a, b, mu, sigma)
		return X
		
		
	def probabilites(self, x, a = None, b = None, mu = 0., sigma=1.):
		'''
			Compute the probabilility of observations
		
		'''
		if a is None:
			a = -np.Inf
		if b is None:
			b = np.Inf
				
		Z = np.sqrt( np.pi/2)*sigma * (erf(b /np.sqrt(2))-erf(a/np.sqrt(2)))
		Z = max(Z, 1e-15)	  # Avoid NaN
		p = np.exp(-(x-mu)**2/2/sigma**2) / Z
		return  p
		
		