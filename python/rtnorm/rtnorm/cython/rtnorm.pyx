import numpy as np
cimport numpy as np
cimport cython



cdef extern from "rtnorm.hpp" namespace "Rtnorm":
	cppclass rtnorm_multi:
		rtnorm_multi(int)
		rtnorm_multi()
		void sample(const int, double* ,const double* ,const double* ,const double* , const double* ) nogil
		
	
cdef class rtnormClass:
	"""
		Class for sampling univariate truncated normal
		
		to generate random variables use sample
	"""
	
	cdef rtnorm_multi* rtnorm_cpp
	
	def __init__(self, Nthreads = -1):
		"""
			Setups the object
		"""
		self.rtnorm_cpp = new rtnorm_multi(Nthreads)
		#TODO: link to a class
		
		
	@cython.boundscheck(False)
	@cython.wraparound(False)			
	def sample(self, np.ndarray[np.double_t, ndim=1] a,
			         np.ndarray[np.double_t, ndim=1] b,
			         np.ndarray[np.double_t, ndim=1] mu,
			         np.ndarray[np.double_t, ndim=1] sigma): 
		
		
		cdef int n = a.shape[0]  # @DuplicatedSignature
		cdef np.ndarray[np.double_t, ndim=1] X = np.zeros((n,))
		with nogil:
			self.rtnorm_cpp.sample(n,
								<double *> &X[0], 
								<double *> &a[0],
								<double *> &b[0], 
								<double *> &mu[0], 
								<double *> &sigma[0])
		#if a.shape > 1
		
		#else
		return X
	

	
		