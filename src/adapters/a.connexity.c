#include <R.h>
#include <Rdefines.h>
#include "sources/connexity.h"

// explore the connectivity of a graph (NIix = neighborhoods indices)
SEXP getConnectedComponents(
	SEXP NIix_
) {
	// extract NIix list vectors in a jagged array
	int n = LENGTH(NIix_);
	int* lengthNIix = (int*)malloc(n*sizeof(int));
	int** NIix = (int**)malloc(n*sizeof(int*));
	for (int i=0; i<n; i++)
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

	// Main call (no R libraries)
	int* connexComps = getConnectedComponents_core(NIix, lengthNIix, n);

	// free memory
	for (int i=0; i<n; i++)
		free(NIix[i]);
	free(NIix);
	free(lengthNIix);

	// transfer result in an R object
	SEXP cc;
	PROTECT(cc = NEW_INTEGER(n));
	int* pcc = INTEGER_POINTER(cc);
	for (int i=0; i<n; i++)
		pcc[i] = connexComps[i];

	// free remaining memory
	free(connexComps);
	UNPROTECT(1);

	return cc;
}
