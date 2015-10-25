//  Example for using rtnorm
//  
//  Copyright (C) 2012 Guillaume Dollé, Vincent Mazet (LSIIT, CNRS/Université de Strasbourg)
//  Licence: GNU General Public License Version 2
//  see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
//
//  Depends: LibGSL
//  OS: Unix based system


#include <iostream>
#include <chrono>
#include <time.h>
#include "rtnorm.hpp"
#define Calloc(n, type) (type *)calloc((size_t)(n),sizeof(type))

using namespace Rtnorm;
int main()
{
  double a = 1;                 // Left bound
  double b = 9;                 // Right bound
  double mu = 2;                // Mean
  double sigma = 3;             // Standard deviation
  double  s;  // Output argument of rtnorm
  int K = 6e5;                  // Number of random variables to generate


  //--- generate and display the random numbers ---
  	rtnorm sampler = rtnorm();
  
    auto wcts = std::chrono::system_clock::now();
    
    for(int k=0; k<K; k++)
  	{
    	s = sampler.sample(a , b, mu, sigma);
  	}
    std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
    
    float  sec_rtnorm = wctduration.count();
  
   	std::random_device rd;
    std::mt19937_64 gen(rd());
	std::normal_distribution<double> N01(0.,1.);
 
  	double r;
  	double a_, b_;
  	for(int k=0; k<5; k++)
  	{
  		a_ = sigma*10*N01(gen);
  		b_ = a_ + sigma*10*fabs(N01(gen));
    	r = sampler.sample(a_,b_, mu, sigma);
    	std::cout << "[ "<< a_ << "," << b_ << "] = ";
    	std::cout<< r <<std::endl;
  	}

	//  NOW REDO WITH THE ONPENMP OBJECT
	
	
	rtnorm_multi samplerOpenmp = rtnorm_multi();
	
	wcts = std::chrono::system_clock::now();
	
	for(int k=0; k<K; k++)
  	{
    	 samplerOpenmp.sample(1, &s, &a , &b, &mu, &sigma);
  	}
  	
  	wctduration = (std::chrono::system_clock::now() - wcts);
    float sec_rtnormM  = wctduration.count();
  	
  	double *x         = Calloc(K,double);
  	double *a_vec     = Calloc(K,double);
  	double *b_vec     = Calloc(K,double);
  	double *mu_vec    = Calloc(K,double);
  	double *sigma_vec = Calloc(K,double);
  	
  	
	for(int k=0; k<K; k++)
  	{
    	a_vec[k]     = a;
    	b_vec[k]     = b;
    	mu_vec[k]    = mu;
    	sigma_vec[k] = sigma;
  	}
    
  	 wcts = std::chrono::system_clock::now();
     
     samplerOpenmp.sample(K, x, a_vec , b_vec, mu_vec, sigma_vec);
  	 
  	 wctduration = (std::chrono::system_clock::now() - wcts);

  
  
  	
  	for(int k = 0; k < 4; k++)
  		std::cout << "x[" << k << "] = " << x[k] << "\n";
  	
  	
    float sec_rtnormM2 = wctduration.count();
  	
	std::cout << "\ntiming result for K = " << K << "\n";
	std::cout<< "               rtnorm       time = "<< sec_rtnorm   << " (sec)" <<std::endl;
	std::cout<< "(poor version) rtnorm_multi time = "<< sec_rtnormM  << " (sec)" <<std::endl;
	std::cout<< "               rtnorm_multi time = "<< sec_rtnormM2 << " (sec)" <<std::endl;

	free(a_vec);
	free(b_vec);
	free(mu_vec);
	free(sigma_vec);
  return 0;
}

