
#include <chrono>
#include <mex.h>
#include <random>
#include "rtnorm.hpp"
#ifdef _OPENMP
#include<omp.h>
#endif

void error_checking(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void value_checking(const int n, const double* a, const double* b, const double* sigma);
/*
-----------------------------------------------------------------------------

MAIN FUNCTION:


-----------------------------------------------------------------------------
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){


//---- defining macros ----//
	#define a_in        prhs[0] // a left limit
	#define b_in        prhs[1] // b right limit
	#define mu_in       prhs[2] // means
	#define sigma_in    prhs[3] // standard devations

	#define X_out       plhs[0] // output variable




//----- error checking -------- //
	error_checking(nlhs, plhs, nrhs, prhs);


//---- linking matlab obj to c obj and checking bad input ---------- //
	
	
	int n_elements =  (int) mxGetNumberOfElements(a_in);
	
	double *a          = mxGetPr(a_in);
	double *b          = mxGetPr(b_in);
	double *mu         = mxGetPr(mu_in);
	double *sigma      = mxGetPr(sigma_in);
	value_checking(n_elements, a, b, mu);
	


	X_out = mxCreateDoubleMatrix(n_elements,
                                 1,
                                 mxREAL);
	double* X = mxGetPr( X_out );


//----sampling ---------- //
	Rtnorm::rtnorm_multi samplerOpenmp = Rtnorm::rtnorm_multi();
	samplerOpenmp.sample(n_elements, X, a , b, mu, sigma);
	//samplerOpenmp.sample(K, x, a_vec , b_vec, mu_vec, sigma_vec);


}




/**Checks if a given array contains a non-sparse, real valued, double matrix.
 *@param tmp pointer to a matlab array.*/
bool mxIsRealDoubleArray(const mxArray* tmp){
  return !mxIsEmpty(tmp) && mxIsDouble(tmp) && !mxIsSparse(tmp) &&
    !mxIsComplex(tmp); 
}

/**Checks if a given array contains a double scalar.
 *@param tmp pointer to a matlab array.*/
bool mxIsRealScalar(const mxArray* tmp){
  return !mxIsEmpty(tmp) && mxIsDouble(tmp) && !mxIsComplex(tmp) &&
    mxGetNumberOfElements(tmp)==1;
}

/**Checks if a given array contains a positive double scalar.
 *@param tmp pointer to a matlab array.*/
bool mxIsRealPosScalar(const mxArray* tmp){
  return mxIsRealScalar(tmp) && mxGetScalar(tmp)>0;
}


void error_checking(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 if( nlhs != 1 ){
    mexErrMsgIdAndTxt("rtnorm:badinput", "Function returns ONE value.");
  }
  if(nrhs != 4 ){
    mexErrMsgIdAndTxt("rtnorm:badinput",
                      "Function requires  4 input(s), %i found.", nrhs);
  }  
  

//check that all four inputs are matrices and extract the sizes
  mwSize sz[3];
  mwSize maxSize=0;
  size_t iMax=0;
  for(size_t i=0; i<4; ++i){
    if( !mxIsRealDoubleArray(prhs[i]) ){
      mexErrMsgIdAndTxt("rtnorm:badinput",
                        "Argument %i should be a matrix", i+1);
    }
    sz[i] = mxGetNumberOfElements(prhs[i]);
    if( sz[i]>maxSize ){
      maxSize=sz[i];
      iMax=i;
    };
  }
  //check that all four inputs have the same size (or are scalar)
  for(size_t i=0; i<4; ++i){
    if(sz[i]!=maxSize ){
      mexErrMsgIdAndTxt("rtnorm:badinput",
                        "Argument %i should have %i elements", i+1, maxSize);
    }
  }

}

void value_checking(const int n,
					const double* a, 
					const double* b,
					const double* sigma)
{

	for(int i = 0; i < n ; i++)
	{
		if( a[i] >= b[i])
			mexErrMsgIdAndTxt("rtnorm:badinput", "a[%i] >= b[%i]",i + 1,i + 1);
	
		if( sigma[i] <= 0)
			mexErrMsgIdAndTxt("rtnorm:badinput", "sigma[%i] <= 0",i  +1);
	}



}