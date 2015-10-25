'''
Created on Oct 24, 2015

@author: jonaswallin
'''

USE_OPENMP = True

import numpy as np
import os
from Cython.Distutils import build_ext
try:
	from setuptools import setup, Extension, find_packages, Command
except ImportError:
	
	try:
		from setuptools.core import setup, Extension, find_packages, Command
	except ImportError:
		from distutils.core import setup, Extension, find_packages, Command




import numpy.distutils.system_info as sysinfo
include_dirs_ = [lib_ for lib_ in sysinfo.default_lib_dirs]

class CleanCommand(Command):
	"""Custom clean command to tidy up the project root."""
	user_options = []
	def initialize_options(self):
		pass
	def finalize_options(self):
		pass
	def run(self):
		os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')



extra_compile_args=["-std=c++11","-O3"]
extra_link_args=["-std=c++11"]
if USE_OPENMP:
	extra_compile_args.append("-fopenmp")
	extra_link_args.append( "-fopenmp")




metadata = dict(
	name='rtnorm',
	cmdclass = {"build_ext": build_ext, 'cleanAll': CleanCommand},
	packages		 = find_packages(),
	package_dir	  = {'rtnorm': 'rtnorm'},
	version		  = '0.1',
	description	  = 'Sampler for truncated normal random variables',
	author		   = 'Jonas Wallin',
	maintainer_email = 'jonas.wallin81@gmail.com',
	url			  = 'https://github.com/JonasWallin/rtnorm',
	author_email	 = 'jonas.wallin81@gmail.com',
	install_requires = ['cython', 'numpy'],
	requires		 = ['numpy (>=1.3.0)','cython'],
	ext_modules	  = [Extension('rtnorm.cython.rtnorm',['rtnorm/cython/rtnorm.pyx',
														'rtnorm/cython/cpp/rtnorm.cpp'],
							library_dirs = include_dirs_,
							include_dirs= [np.get_include(),'rtnorm/cython/cpp/'],
						    extra_compile_args = extra_compile_args,
						    extra_link_args    = extra_link_args,
							language="c++")]
	)
setup(**metadata)