#ifndef SYNCLUST_ALGEBRA_H
#define SYNCLUST_ALGEBRA_H

// small useful function to transform a matrix as given by R 
// into a easier-to-handle one.
double* transpose(
	double* M, 
	int nrow, 
	int ncol
);

// auxiliary to compute euclidian norm
double norm2(
	double* v, 
	int length
);

// auxiliary to compute euclidian distance
double distance2(
	double* v1, 
	double* v2, 
	int length
);

#endif
