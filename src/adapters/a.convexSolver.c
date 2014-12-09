#include <R.h>
#include <Rdefines.h>
#include "sources/convexSolver.h"
#include "sources/utils/algebra.h"

// compute estimated ("repaired", "smoothed"...) variations from rows of M
// NOTE: geographic coordinates dropped here, since they are unused
SEXP getVarsWithConvexOptim(
	SEXP M_, 
	SEXP NIix_, 
	SEXP alpha_, 
	SEXP h_, 
	SEXP epsilon_, 
	SEXP maxiter_, 
	SEXP symmNeighbs_, 
	SEXP trace_
) {
	// get parameters
	double alpha = NUMERIC_VALUE(alpha_);
	double h = NUMERIC_VALUE(h_);
	double epsilon = NUMERIC_VALUE(epsilon_);
	int maxiter = INTEGER_VALUE(maxiter_);
	int symmNeighbs = LOGICAL_VALUE(symmNeighbs_);
	int trace = LOGICAL_VALUE(trace_);

	// extract infos from M and get associate pointer
	SEXP dim = getAttrib(M_, R_DimSymbol);
	int nrow = INTEGER(dim)[0];
	int ncol = INTEGER(dim)[1];
	// M is always given by columns: easier to process in rows
	double* pM = transpose(REAL(M_), nrow, ncol);

	// extract NIix list vectors in a jagged array
	int* lengthNIix = (int*)malloc(nrow*sizeof(int));
	int** NIix = (int**)malloc(nrow*sizeof(int*));
	for (int i=0; i<nrow; i++)
	{
		lengthNIix[i] = LENGTH(VECTOR_ELT(NIix_,i));
		SEXP tmp;
		PROTECT(tmp = AS_INTEGER(VECTOR_ELT(NIix_,i)));
		NIix[i] = (int*)malloc(lengthNIix[i]*sizeof(int));
		for (int j=0; j<lengthNIix[i]; j++)
			NIix[i][j] = INTEGER(tmp)[j];
		UNPROTECT(1);
		// WARNING: R indices start at 1,
		// so we must lower every index right now to avoid future bug
		for (int j=0; j<lengthNIix[i]; j++)
			NIix[i][j]--;
	}

	// Main call to core algorithm
	Parameters params = getVarsWithConvexOptim_core(
		pM, lengthNIix, NIix, nrow, ncol, alpha, h, epsilon, maxiter, symmNeighbs, trace);

	// free neighborhoods parameters arrays
	free(lengthNIix);
	for (int i=0; i<nrow; i++)
		free(NIix[i]);
	free(NIix);

	// copy matrix F into pF for output to R (1D matrices)
	SEXP f;
	PROTECT(f = allocMatrix(REALSXP, nrow, ncol));
	double* pF = REAL(f);
	for (int i=0; i<nrow; i++)
	{
		for (int j=0; j<ncol; j++)
			pF[i+nrow*j] = params.f[i][j];
	}
	// copy theta into pTheta for output to R
	SEXP theta;
	PROTECT(theta = allocVector(REALSXP, nrow));
	double* pTheta = REAL(theta);
	for (int i=0; i<nrow; i++)
		pTheta[i] = params.theta[i];

	// free params.f and params.theta
	free(params.theta);
	for (int i=0; i<nrow; i++)
		free(params.f[i]);
	free(params.f);

	// build return list with f and theta
	SEXP listParams, listNames;
	PROTECT(listParams = allocVector(VECSXP, 2));
	char* lnames[2] = {"f", "theta"}; //lists labels
	PROTECT(listNames = allocVector(STRSXP,2));
	for (int i=0; i<2; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));
	setAttrib(listParams, R_NamesSymbol, listNames);
	SET_VECTOR_ELT(listParams, 0, f);
	SET_VECTOR_ELT(listParams, 1, theta);

	UNPROTECT(4);
	return listParams;
}
