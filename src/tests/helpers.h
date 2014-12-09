#ifndef SYNCLUST_HELPERS_H
#define SYNCLUST_HELPERS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "sources/utils/boolean.h"

//auxiliary to check vectors equality
int checkEqualV(
	double* v1, 
	double* v2, 
	int n
);

// auxiliary to count distinct values in an integer array
int countDistinctValues(
	int* v, 
	int n
);

// check if clusters proportions match a given fraction (tolerance 'tol')
int checkClustersProportions(
	int* clusters, 
	int n, 
	int clustCount, 
	double tol
);

// basic unit tests macros
#define ASSERT_TRUE(b) \
do { \
	if (!(b)) { \
		fprintf(stdout, "Error in file %s at line %i\n", __FILE__, __LINE__); \
		return; \
	} \
} while (0)

#define ASSERT_FALSE(b) \
do { \
	if (b) { \
		fprintf(stdout, "Error in file %s at line %i\n", __FILE__, __LINE__); \
		return; \
	} \
} while (0)

#endif
