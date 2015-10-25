# rtnorm
Sampling truncated univ normal distribution.


This repos contains python and matlab package for simulation of 
truncated normal random variables based upon the codes:
[Vincent Mazet](http://miv.u-strasbg.fr/mazet/rtnorm/), and it built from 
the article [Fast simulation of trunacted Gaussian distributions](http://link.springer.com/article/10.1007%2Fs11222-009-9168-1) of [Niclos Chopin](https://sites.google.com/site/nicolaschopinstatistician/software).



## Python
For python the requirements are:

* scipy (will be removed)
* NumPy
* Cython
* openmp (option can be turned off in setupy.py)



		pip -e 'git+https://git@github.com/JonasWallin/rtnorm.git#egg=rtnorm&subdirectory=python/rtnorm' 

To use the sampler:

	import rtnorm
	sampler = rtnorm.rtnorm()
	X = sampler(a = 1, b = 2, mu = 3, sigma = 4)
	
Or with vectors:

	import rtnorm
	import numpy as np
	X = sampler( a = np.zeros((4,)), b = np.random.rand(4) )



## Matlab

For matlab: 

	make 


should install the make mex files link to the directory with addpath.
Run matlab code with:

	x = rtnorm(a, b, mu, sigma);

where a,b,mu, and sigma are vectors

###TODO:

	1. setup test case  in python (done)
	2. setup speed test in python (done)
	3. modify the code so it sutiable for vector versions (done)
	4. Write the probability so we can test more advanced 
	5. Matlab version (done)
	7. test on the amzaon cluster (done)