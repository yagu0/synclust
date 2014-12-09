#include <R.h>
#include <Rdefines.h>
#include "sources/dijkstra.h"

// Dijkstra from index start : return vector of distances to every other vertex
// NOTE [space VS perf]: passing neighborhoods infos would be enough, but
//					   require extra computation to translate R list argument
SEXP dijkstra(
	SEXP distances_, 
	SEXP start_
) {
	// get arguments
	SEXP dim = getAttrib(distances_, R_DimSymbol);
	int n = INTEGER(dim)[0];
	double* pDistsIn = REAL(distances_);
	int start = INTEGER_VALUE(start_) - 1; // R indices start at 1

	// Main call to core algorithm
	double* distances = dijkstra_core(pDistsIn, start, n);

	// allocate vector output and obtain R vector object
	SEXP distsOut;
	PROTECT(distsOut = allocVector(REALSXP, n));
	double* pDistsOut = NUMERIC_POINTER(distsOut);
	for (int i=0; i<n; i++)
		pDistsOut[i] = distances[i];

	free(distances);
	UNPROTECT(1);

	return distsOut;
}
