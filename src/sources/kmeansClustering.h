#ifndef SYNCLUST_KMEANSCLUSTERING_H
#define SYNCLUST_KMEANSCLUSTERING_H

#include <cgds/Vector.h>

// auxiliary function to obtain a random sample of 1..n with K elements
void sample(
	int* centers, 
	int n, 
	int K
);

// auxiliary function to compare two sets of centers
int unequalCenters(
	int* ctrs1, 
	int* ctrs2, 
	int n, 
	int K
);

// assign a vector (represented by its distances to others, as distances[index,])
// to a cluster, represented by its center index ==> output is integer in 0..K-1
int assignCluster(
	int index, 
	double* distances, 
	int* centers, 
	int n, 
	int K
);

// k-means based on a distance matrix (nstart=10, maxiter=100)
int* kmeansWithDistances_core(
	double* pDistances, 
	int n, 
	int K, 
	int nstart, 
	int maxiter
);

#endif
