#include "sources/utils/algebra.h"
#include <stdlib.h>
#include <math.h>

// small useful function to transform a matrix as given by R 
// into a easier-to-handle one.
double* transpose(double* M, int nrow, int ncol)
{
	double* Mtr = (double*)malloc(nrow*ncol*sizeof(double));
	for (int i=0; i<nrow; i++)
	{
		for (int j=0; j<ncol; j++)
			Mtr[i*ncol+j] = M[i+nrow*j];
	}
	return Mtr;
}

// auxiliary to compute euclidian norm
double norm2(double* v, int length)
{
	double norm = 0.0;
	for (int j=0; j<length; j++)
		norm += v[j]*v[j];
	return sqrt(norm);
}

// auxiliary to compute euclidian distance
double distance2(double* v1, double* v2, int length)
{
	double distance = 0.0, diff;
	for (int j=0; j<length; j++)
	{
		diff = v1[j]-v2[j];
		distance += diff*diff;
	}
	return sqrt(distance);
}
