#ifndef SYNCLUST_CONVEXSOLVER_H
#define SYNCLUST_CONVEXSOLVER_H

#include "sources/utils/boolean.h"

// auxiliary to compute euclidian norm
double norm2(
	double* v, 
	int length
);

// auxiliary to compute euclidian distance
double distance2(
	double* f1, 
	double* f2, 
	int length
);

// auxiliary to compute log-likelihood + penalty
double computeLogLikelihood(
	double** f, 
	double* theta, 
	double** Zst, 
	double*** phi, 
	int* lengthNIix, 
	int** NIix, 
	double alpha, 
	int nrow, 
	int ncol
);

// structure to return parameters (theta, f) [and others if needed later]
typedef struct Parameters {
    double** f;
    double* theta;
} Parameters;

// compute estimated ("repaired", "smoothed"...) variations from rows of M
// NOTE: geographic coordinates dropped here, since they are unused
Parameters getVarsWithConvexOptim_core(
	double* pM, 
	int* lengthNIix, 
	int** NIix, 
	int nrow, 
	int ncol, 
	double alpha, 
	double h, 
	double epsilon, 
	int maxiter, 
	bool symmNeighbs, 
	bool trace
);

#endif
