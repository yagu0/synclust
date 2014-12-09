#include <R.h>
#include <Rdefines.h>
#include "sources/kmeansClustering.h"
#include <cgds/Vector.h>

// k-means based on a distance matrix (nstart=10, maxiter=100)
SEXP kmeansWithDistances(
	SEXP distances_, 
	SEXP K_, 
	SEXP nstart_, 
	SEXP maxiter_
) {
	// get scalar arguments
	int K = INTEGER_VALUE(K_);
	int nstart = NUMERIC_VALUE(nstart_);
	int maxiter = INTEGER_VALUE(maxiter_);

	// extract infos from M and get associate pointer
	SEXP dim = getAttrib(distances_, R_DimSymbol);
	int n = INTEGER(dim)[0];
	double* pDistances = REAL(distances_);

	// Main call to core algorithm
	int* clusters = kmeansWithDistances_core(pDistances, n, K, nstart, maxiter);

	// allocations and recopies to R vector object
	SEXP bestClusts;
	PROTECT(bestClusts = allocVector(INTSXP, n));
	int* pBestClusts = INTEGER(bestClusts);
	for (int i=0; i<n; i++)
		pBestClusts[i] = clusters[i] + 1; // add 1 to start labels at 1
	free(clusters);

	// and return clusters
	UNPROTECT(1);
	return bestClusts;
}
